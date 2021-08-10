/*
    Copyright (C) 2017 Daniel Schultz
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

/* #define EXPTEST (*exp == 0x6000002000004LL) */
#define EXPTEST 0
#define DEBUGTEST 0

typedef unsigned short hind_t;
#define HIND_FORMAT "%d"

/* TMP_ALLOC will use alloca for allocations of 8192 or less, which causes problems */
/* Allocate chain in 16KB blocks */
static const int chain_block_size = 16384/sizeof(mpoly_heap_t);

/*
   Set A to B1*B2*---*Bm + Bm+1*Bm+2*---*Bn + ... using Johnson's heap
   method. The function reallocates its output and returns the length
   of the product. This version of the function assumes the exponent
   vectors take N words.

   mul_johnson.c optimizes for medium sized coefficients (the coeffs
   fit in a single machine word but their products do not) by keeping
   the intermediate results in local variables until we're ready to
   store them in memory.  This routine does not currently implement
   such an optimization, though it could.

   This version of the addmul routine is designed for threading.

   We only want to output a block of data at a time.
*/

/* per-block data structure accessed by the block threads */

struct _fmpz_mpoly_addmul_state {
   const fmpz_mpoly_struct * Blist;
   slong Blength;
   slong num_malloced_memory_blocks;
   void ** malloced_memory_blocks;
   slong next_loc;
   slong heap_len;
   slong heap_size;
   slong heap_block_size;
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** chain_list;
   ulong * exps;
   ulong ** exp_list;
   slong exp_next;
   slong chain_next;
   slong chain_size;
   hind_t * hind;
   slong hind_len;
};

/* per-block data structure that will accessed mainly by the merge thread,
 * and is therefore kept separate from the state structure, which will
 * be accessed mainly by the block threads.
 */

struct _fmpz_mpoly_addmul_multi_control
{
    fmpz * coeffs;
    ulong * exps;
    struct _fmpz_mpoly_addmul_state * state;

    /* everything_generated - a boolean indicating that all of the items from this term have been generated
     * total_generated - a count of the number of items generated
     * total_transferred - a count of the number of items transferred to the merge heap
     *         (note that the merge heap only contains indices into the blocks, so "transferred" is a misnomer)
     * total_output - a count of the number of items output to the final polynomial
     *
     * process_thread is NULL if no thread is generating for this term
     */
    slong everything_generated;
    slong total_generated;
    slong total_transferred;
    slong total_output;
    const struct _fmpz_mpoly_addmul_multi_worker * process_thread;
};

/* heap data structure */

struct _fmpz_mpoly_addmul_multi_heap
{
    short unsigned int control;
    short unsigned int index;
};

/* master data structure for an entire addmul_multi computation */

struct _fmpz_mpoly_addmul_multi_master
{
    fmpz_mpoly_struct * A;
    slong k;
    const fmpz_mpoly_struct * Blist;
    const slong * Blengths;
    slong numterms;
    slong Btotallen;
    flint_bitcnt_t bits;
    slong N;
    const ulong * cmpmask;
    const fmpz_mpoly_ctx_struct * ctx;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    struct _fmpz_mpoly_addmul_multi_control * control;
    struct _fmpz_mpoly_addmul_multi_heap * heap;
    slong heaplen;
    const struct _fmpz_mpoly_addmul_multi_worker * merge_thread;
};

static int blocksize = 1024;
static int numblocks = 4;

/* argument passed to a worker thread
 *
 * We can only pass one argument to each thread.
 * So all of the structures have identical content, but we use the address of the
 * worker structure to distinguish between the various threads.
 */

struct _fmpz_mpoly_addmul_multi_worker
{
    struct _fmpz_mpoly_addmul_multi_master * master;
};

/* Initialize state for processing a single product */

void _fmpz_mpoly_addmul_multi_init_state(
    const fmpz_mpoly_struct * Blist,
    const slong Blength,
    struct _fmpz_mpoly_addmul_state * state,
    struct _fmpz_mpoly_addmul_multi_master * master)
{
   slong i, j;
   mpoly_heap_t * x;

   state->Blist = Blist;
   state->Blength = Blength;

   /* First polynomial should be the largest; it's left out of this calculation */
   state->hind_len = 1;

   /* XXX - FLINT_ASSERTs are only tested if configured with --enable-assert */
   FLINT_ASSERT(Blist[0].length < (UWORD(1) << (8*sizeof(hind_t) - 1)));

   for (j = 1; j < Blength; j++)
   {
       FLINT_ASSERT(Blist[j].length <= Blist[0].length);
       state->hind_len *= Blist[j].length;
   }

   state->next_loc = state->hind_len + 4;   /* something bigger than heap can ever be */

   state->heap_block_size = 16384/sizeof(ulong)/master->N;
   state->heap_size = state->heap_block_size;
   state->chain_size = chain_block_size;

   state->heap = (mpoly_heap_s *) flint_malloc(state->heap_size*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   state->chain = (mpoly_heap_t *) flint_malloc(state->chain_size*sizeof(mpoly_heap_t));
   state->chain_list = (mpoly_heap_t **) flint_malloc(state->chain_size*sizeof(mpoly_heap_t **));
   for (i = 0; i < state->chain_size; i++)
      state->chain_list[i] = state->chain + i;

   /* allocate space for exponent vectors of N words */
   state->exps = (ulong *) flint_malloc(state->heap_size * master->N * sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   state->exp_list = (ulong **) flint_malloc(state->heap_size*sizeof(ulong *));
   for (i = 0; i < state->heap_size; i++)
      state->exp_list[i] = state->exps + i*master->N;

   state->malloced_memory_blocks = (void **) flint_malloc(2 * sizeof(void *));
   state->malloced_memory_blocks[0] = state->chain;
   state->malloced_memory_blocks[1] = state->exps;
   state->num_malloced_memory_blocks = 2;

   /* space for heap indices */
   state->hind = (hind_t *) flint_malloc(state->hind_len*sizeof(hind_t));
   for (i = 0; i < state->hind_len; i++)
       state->hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   state->heap_len = 1; /* heap starts empty, and its zero index is unused, so heap_len = 1 */
   state->exp_next = 0;
   state->chain_next = 0;

   /* put (0, exp2[0] + ... + expn[0]) on heap */
   x = state->chain_list[state->chain_next ++];
   x->i = 0;
   x->next = NULL;

   mpoly_monomial_zero(state->exp_list[state->exp_next], master->N);
   for (i = 0; i < Blength; i++)
       if (master->bits <= FLINT_BITS)
           mpoly_monomial_add(state->exp_list[state->exp_next], state->exp_list[state->exp_next], Blist[i].exps, master->N);
       else
           mpoly_monomial_add_mp(state->exp_list[state->exp_next], state->exp_list[state->exp_next], Blist[i].exps, master->N);

   state->hind[0] = 2*1 + 0;

   if (!_mpoly_heap_insert(state->heap, state->exp_list[state->exp_next++], x,
                           &state->next_loc, &state->heap_len, master->N, master->cmpmask))
      state->exp_next--;

}

slong _fmpz_mpoly_addmul_multi_process_block(
    fmpz * coeffs,
    ulong * exps,
    slong block_size,
    struct _fmpz_mpoly_addmul_state * state,
    struct _fmpz_mpoly_addmul_multi_master * master)
{
   slong i, j, k, l;
   mpoly_heap_t * x;
   ulong * exp;
   ulong * Q;
   slong Q_len;
   slong Q_size;
   ulong multiindex;
   int first;
   ulong offset, offset2, offset3;
   ulong candidate;
   ulong partial_multiindex;
   fmpz_t tmp_coeff;

   fmpz_init(tmp_coeff);

   /* space for temporary storage of pointers to heap nodes */
   Q_size = 16384/sizeof(ulong);
   Q = (ulong *) flint_malloc(Q_size*sizeof(ulong));
   Q_len = 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while ((state->heap_len > 1) && (k+1 < block_size))
   {
      /* get pointer to exponent field of heap top */
      exp = state->heap[1].exp;

      /* realloc output poly ready for next product term */
      k++;

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* while heap nonempty and contains chain with current output exponent */
      while (state->heap_len > 1 && mpoly_monomial_equal(state->heap[1].exp, exp, master->N))
      {
         /* pop chain from heap and set exponent field to be reused */
         state->exp_list[--state->exp_next] = state->heap[1].exp;

         x = _mpoly_heap_pop(state->heap, &state->heap_len, master->N, master->cmpmask);

         if (first)
         {
            fmpz_zero(coeffs + k);
            
            /* set output monomial */
            mpoly_monomial_set(exps + k*master->N, exp, master->N);

            first = 0; 
         }

         /* for each node in this chain */
         do
         {
            multiindex = x->i;

            if (EXPTEST) {
                fprintf(stderr, "Processing %ld with exp %lx\n", multiindex, *exp);
            }
            /* addmul product of input poly coeffs */

            fmpz_set(tmp_coeff, state->Blist[0].coeffs + (state->hind[multiindex] >> 1) - 1);
            partial_multiindex = multiindex;
            for (i=1; i<state->Blength; i++)
            {
                fmpz_mul(tmp_coeff, tmp_coeff, state->Blist[i].coeffs + (partial_multiindex % state->Blist[i].length));
                partial_multiindex /= state->Blist[i].length;
            }
            fmpz_add(coeffs + k, coeffs + k, tmp_coeff);

            if (EXPTEST) {
                fprintf(stderr, "adding tmp_coeff = ");
                fmpz_fprint(stderr, tmp_coeff);
                fprintf(stderr, " to get A->coeffs[%ld]=", k);
                fmpz_fprint(stderr, coeffs + k);
                fprintf(stderr, "\n");
            }

            /* take node out of heap and put into store */
            if (Q_len == Q_size)
            {
                Q_size += 16384/sizeof(ulong);
                Q = (ulong *) flint_realloc(Q, Q_size*sizeof(ulong));
            }
            state->hind[multiindex] |= WORD(1);
            Q[Q_len++] = multiindex;

            state->chain_list[-- state->chain_next] = x;
         }
         while ((x = x->next) != NULL);
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         multiindex = Q[--Q_len];

         offset = 0;
         for (i=0; i<state->Blength; i++)
         {
             candidate = multiindex + offset;
             if (EXPTEST) fprintf(stderr, "Considering candidate %ld\n", candidate);
             if ((candidate < state->hind_len) && (state->hind[candidate] < 2*state->Blist[0].length) && (state->hind[candidate] & 1))
             {
                 offset2 = 1;
                 for (j=1; j<state->Blength; j++)
                 {
                     offset3 = offset2 * state->Blist[j].length;
                     if (((candidate % offset3) / offset2 != 0) && (state->hind[candidate - offset2] < state->hind[candidate] + 2))
                         break;
                     offset2 *= state->Blist[j].length;
                 }
                 if (j == state->Blength)
                 {
                     if (state->heap_len == state->heap_size)
                     {
                         state->heap_size += state->heap_block_size;
                         state->heap = flint_realloc(state->heap, state->heap_size*sizeof(mpoly_heap_s));
                         state->exps = flint_malloc(state->heap_block_size*master->N*sizeof(ulong));
                         state->exp_list = flint_realloc(state->exp_list, state->heap_size*sizeof(ulong *));
                         for (l = 0; l < state->heap_block_size; l++)
                             state->exp_list[state->heap_len + l] = state->exps + l*master->N;

                         state->num_malloced_memory_blocks ++;
                         state->malloced_memory_blocks = (void **) flint_realloc(state->malloced_memory_blocks,
                                                                                  state->num_malloced_memory_blocks * sizeof(void *));
                         state->malloced_memory_blocks[state->num_malloced_memory_blocks - 1] = state->exps;
                     }

                     if (state->chain_next == state->chain_size)
                     {
                         state->chain_size += chain_block_size;
                         state->chain = (mpoly_heap_t *) flint_malloc(chain_block_size*sizeof(mpoly_heap_t));
                         state->chain_list = (mpoly_heap_t **) flint_realloc(state->chain_list, state->chain_size*sizeof(mpoly_heap_t *));
                         for (l = 0; l < chain_block_size; l++)
                             state->chain_list[state->chain_next + l] = state->chain + l;

                         state->num_malloced_memory_blocks ++;
                         state->malloced_memory_blocks = (void **) flint_realloc(state->malloced_memory_blocks,
                                                                                  state->num_malloced_memory_blocks * sizeof(void *));
                         state->malloced_memory_blocks[state->num_malloced_memory_blocks - 1] = state->chain;
                     }

                     x = state->chain_list[state->chain_next ++];
                     x->i = candidate;
                     x->next = NULL;

                     mpoly_monomial_set(state->exp_list[state->exp_next], state->Blist[0].exps + master->N*(state->hind[candidate] >> 1), master->N);
                     state->hind[candidate] ++;

                     partial_multiindex = candidate;
                     for (l = 1; l < state->Blength; l++)
                     {
                         if (master->bits <= FLINT_BITS)
                             mpoly_monomial_add(state->exp_list[state->exp_next], state->exp_list[state->exp_next],
                                                state->Blist[l].exps + master->N*(partial_multiindex % state->Blist[l].length), master->N);
                         else
                             mpoly_monomial_add_mp(state->exp_list[state->exp_next], state->exp_list[state->exp_next],
                                                   state->Blist[l].exps + master->N*(partial_multiindex % state->Blist[l].length), master->N);
                         partial_multiindex /= state->Blist[l].length;
                     }

                     if (EXPTEST) {
	                 fprintf(stderr, "Adding %ld because of %ld with hind[%ld] " HIND_FORMAT "\n",
                                 candidate, multiindex, candidate, state->hind[candidate]);
                     }
                     if (!_mpoly_heap_insert(state->heap, state->exp_list[state->exp_next++], x,
                                             &state->next_loc, &state->heap_len, master->N, master->cmpmask))
                         state->exp_next--;
                 }
             }

             if (offset == 0)
                 offset = 1;
             else
                 offset *= state->Blist[i].length;
         }
      }

      if (fmpz_is_zero(coeffs + k))
         k--;
   }

   k++;

   fmpz_clear(tmp_coeff);
   flint_free(Q);

   return k;
}

/* For each term in the sum, we form a struct _fmpz_mpoly_addmul_state
 * and two 4-blocks (one for coeffs and one for exps).  We also keep
 * the index of the next term to remove from the blocks and the index
 * of the last valid term.  I want to keep that data separate from
 * _fmpz_mpoly_addmul_state to keep it on a different cache line
 * (64 bytes most commonly) since different threads will access
 * the two different structures.
 *
 *
 * We also need a heap structure that feeds from those blocks.
 *
 * heap structure: pointer to exponent vector in exp block
 *                 pointer to beginning of exp block
 *                 index into exp/coeff blocks
 *                 pointer to struct with pointers to exp/coeff blocks
 *                 index to struct with pointers to exp/coeff blocks
 *                 min 8 bits index to struct, since some of our addmuls have 130 input terms
 *               * short index to struct and short index into blocks
 *
 * struct control
 *    fmpz * coeffs  (4096 words)
 *    ulong * exps   (N*4096 words)
 *    addmul_state *
 *    
 *
 * We wait until all of the 1-blocks have been created before
 * initializing the heap structure.  Then we have a subroutine
 * that pulls from the heap structure until it stalls and writes
 * an output polynomial.
 *
 * We initialize all of the addmul_state structures.  We run
 * process_block once for each term.  Then we initialize the
 * heap structure and run the subroutine that feeds from it.
 *
 * When a process_block finishes, we atomically adjust the
 * last valid term, which may allow the output routine to keep
 * running.  Otherwise, when the output routine stalls, we
 * run an extra process_block.  If we fill up all of the
 * 4-blocks and can't run any more process_blocks until
 * the output routine finishes, we print a warning message.
 */

void _fmpz_mpoly_addmul_multi_merge(
    struct _fmpz_mpoly_addmul_multi_master * master
)
{
    const slong N = master->N;
    const ulong * cmpmask = master->cmpmask;

    slong i, j;
    int need_to_block = 0;
    int run_until_heap_consumed = 1;
    int first;
    struct _fmpz_mpoly_addmul_multi_heap * heap = master->heap;
    struct _fmpz_mpoly_addmul_multi_control * control = master->control;

    /* start by filling in the end of the heap with data from the blocks */

    for (i=0; i < master->numterms; i++)
    {
        if (heap[master->heaplen/2 + i].index == (unsigned short int ) -1)
        {
            if (control[i].total_generated > control[i].total_transferred)
            {
                heap[master->heaplen/2 + i].index = control[heap[1].control].total_transferred % (blocksize * numblocks);
                control[heap[1].control].total_transferred ++;
            }
        }
        /* If there might be more generated data from any term, we can't run_until_heap_consumed */
        if (! control[i].everything_generated)
            run_until_heap_consumed = 0;
    }

    while (! need_to_block && (heap[1].index != (unsigned short int) -1))
    {

        /* since the process_block routine combines all identical exponents in each block, */
        /* we can be sure that all identical exponents are present in the heap so long */
        /* as each block has been loaded */

        master->k ++;
        fmpz_mpoly_fit_length(master->A, master->k + 1, master->ctx);

        mpoly_monomial_set(master->A->exps + master->k*master->N, control[heap[1].control].exps + master->N*heap[1].index, master->N);
        first = 1;

        /* XXX need to keep going if need_to_block; I put this here for the 1-block case (heap[1] hasn't been updated) */
        while ((heap[1].index != (unsigned short int) -1) && mpoly_monomial_equal(master->A->exps + master->k*N, control[heap[1].control].exps + master->N*heap[1].index, master->N))
        {
            if (first)
                fmpz_swap(master->A->coeffs + master->k, control[heap[1].control].coeffs + heap[1].index);
            else
                fmpz_add(master->A->coeffs + master->k, master->A->coeffs + master->k, control[heap[1].control].coeffs + heap[1].index);

            first = 0;

            control[heap[1].control].total_output ++;

            /* pop lead item from heap and propagate items through the heap */
            i = 1;
            j = 2;
            while (j < master->heaplen)
            {
                /* both j and j+1 have hit their ends */
                if ((heap[j].index == (unsigned short int) -1) && (heap[j+1].index == (unsigned short int) -1))
                {
                    heap[i].index = -1;
                    break;
                }
                /* take the larger of j and j+1, or the one that hasn't hit its end */
                if ((heap[j].index == (unsigned short int) -1)
                    || ((heap[j+1].index != (unsigned short int) -1)
                        && !mpoly_monomial_gt(control[heap[j].control].exps + master->N*heap[j].index,
                                              control[heap[j+1].control].exps + master->N*heap[j+1].index, N, cmpmask)))
                    j ++;
                heap[i] = heap[j];
                i = j;
                j = HEAP_LEFT(j);
            }

            if ((i >= master->heaplen/2) && (heap[i].index != (unsigned short int) -1))
            {
                if (control[heap[1].control].total_transferred == control[heap[1].control].total_generated)
                {
                    heap[i].index = -1;
                    if (! run_until_heap_consumed)
                        need_to_block = 1;
                }
                else
                {
                    heap[i].index ++;
                    heap[i].index %= (blocksize * numblocks);
                    control[heap[1].control].total_transferred ++;
                }
            }
        }
    }
}

void _fmpz_mpoly_addmul_multi_merge_init(
    struct _fmpz_mpoly_addmul_multi_master * master
)
{
    slong i,j,k,l;
    struct _fmpz_mpoly_addmul_multi_control * control = master->control;
    struct _fmpz_mpoly_addmul_multi_heap * heap;
    slong heaplen;

    for (k = 0; ; k++) /* round up numterms to a power of 2 */
    {
        if (master->numterms <= (1<<k))
            break;
    }

    /* (1<<k) is the number of input slots we need at the end of the heap */
    /* (2<<k) is the required size of the heap */

    heaplen = 2<<k;
    heap = (struct _fmpz_mpoly_addmul_multi_heap *) flint_calloc(heaplen, sizeof(struct _fmpz_mpoly_addmul_multi_heap));

    for (i=0; i<heaplen; i++)
    {
        if (i < heaplen/2)
        {
            heap[i].control = 0;
            heap[i].index = -1;
        }
        else if (i - heaplen/2 < master->numterms)
        {
            heap[i].control = i - heaplen/2;
            heap[i].index = 0;
            control[i - heaplen/2].total_transferred ++;
        }
        else
        {
            heap[i].control = 0;
            heap[i].index = -1;
        }
    }

    for (l=0; l<k; l++)
    {
        for (i=1; i<heaplen/2; i++)
        {
            j = HEAP_LEFT(i);
            if ((heap[j].index == (unsigned short int) -1)
                || ((heap[j+1].index != (unsigned short int) -1)
                    && !mpoly_monomial_gt(master->control[heap[j].control].exps + master->N*heap[j].index,
                                          master->control[heap[j+1].control].exps + master->N*heap[j+1].index, master->N, master->cmpmask)))
                j ++;
            heap[i] = heap[j];

            if ((j >= heaplen/2) && (j - heaplen/2 < master->numterms) && (heap[j].index != (unsigned short int) -1))
            {
                if (control[j - heaplen/2].total_transferred == control[j - heaplen/2].total_generated)
                {
                    /* Assume that the input blocks are larger than the log of the heap size,
                     * so if we run out of data it's because there is no more, not because
                     * we need to fetch another block.
                     */
                    control[j - heaplen/2].everything_generated = UWORD(1);
                    heap[j].index = -1;
                }
                else
                {
                    control[j - heaplen/2].total_transferred ++;
                    heap[j].index ++;
                }
            }
        }
    }

    FLINT_ASSERT(heap[0].index == (unsigned short int) -1);
    FLINT_ASSERT(heap[1].index == 0);

    master->heap = heap;
    master->heaplen = heaplen;
}

static void _fmpz_mpoly_addmul_multi_threaded_worker(void * varg)
{
    struct _fmpz_mpoly_addmul_multi_worker * worker = (struct _fmpz_mpoly_addmul_multi_worker *) varg;
    struct _fmpz_mpoly_addmul_multi_master * master = worker->master;
    slong i, j, k;
    slong offset;
    slong able_to_merge;
    slong all_done;
    slong least_valid_data;
    slong least_valid_data_block;

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(& master->mutex);
#endif
    while (1)
    {
        able_to_merge = 1;
        all_done = 1;
        /* We want to find the term with the least amount of data generated but not output,
         * but we don't want to pick any term if there isn't at least one free block available,
         * so start with (numblocks - 1)*blocksize and look for terms with less data than that.
         */
        least_valid_data = (numblocks - 1)*blocksize;
        least_valid_data_block = -1;
        for (i = 0; i < master->numterms; i ++)
        {
            FLINT_ASSERT(master->control[i].total_output <= master->control[i].total_transferred);
            FLINT_ASSERT(master->control[i].total_transferred <= master->control[i].total_generated);
            if ((! master->control[i].everything_generated) || (master->control[i].total_output != master->control[i].total_generated))
                all_done = 0;
            if ((! master->control[i].everything_generated) && (master->control[i].total_generated == master->control[i].total_transferred))
                able_to_merge = 0;
            if ((! master->control[i].everything_generated) && (! master->control[i].process_thread))
            {
                slong valid_data = (master->control[i].total_generated - master->control[i].total_output);
                if (valid_data < least_valid_data)
                {
                    least_valid_data = valid_data;
                    least_valid_data_block = i;
                }
            }
        }

        all_done = all_done && (master->heap != NULL) && (master->heap[1].index == (unsigned short int) -1);

        if (!all_done && able_to_merge && (! master->merge_thread))
        {
            /* start a merge running; */
            if (master->heap == NULL)
                _fmpz_mpoly_addmul_multi_merge_init(master);
            master->merge_thread = worker;
            if (DEBUGTEST) fprintf(stderr, "threaded_worker %p merging\n", varg);
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(& master->mutex);
#endif
            _fmpz_mpoly_addmul_multi_merge(master);
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(& master->mutex);
#endif
            master->merge_thread = NULL;
        }
        else if (!all_done && (least_valid_data_block >= 0))
        {
            i = least_valid_data_block;
            master->control[i].process_thread = worker;
            if (DEBUGTEST) fprintf(stderr, "threaded_worker %p block %ld\n", varg, i);
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(& master->mutex);
#endif
            if (master->control[i].state == NULL)
            {
                master->control[i].state = (struct _fmpz_mpoly_addmul_state *) flint_malloc(sizeof(struct _fmpz_mpoly_addmul_state));
                j = 0;
                for (k = 0; k < i; k ++)
                    j += master->Blengths[i];
                _fmpz_mpoly_addmul_multi_init_state(master->Blist+j, master->Blengths[i], master->control[i].state, master);
                master->control[i].exps = (ulong *) flint_malloc(4 * blocksize * master->N * sizeof(ulong));
                master->control[i].coeffs = (fmpz *) flint_malloc(4 * blocksize * sizeof(fmpz));
                for (k = 0; k < 4 * blocksize; k ++)
                    fmpz_init(master->control[i].coeffs + k);
            }
            offset = master->control[i].total_generated % (blocksize * numblocks);
            k = _fmpz_mpoly_addmul_multi_process_block(master->control[i].coeffs + offset,
                                                       master->control[i].exps + master->N * offset,
                                                       blocksize, master->control[i].state, master); /*  */
            if (DEBUGTEST) fprintf(stderr, "threaded_worker %p block %ld returns %ld\n", varg, i, k);
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(& master->mutex);
#endif
            if (k)
                master->control[i].total_generated += k;
            else
                master->control[i].everything_generated = UWORD(1);
            master->control[i].process_thread = NULL;
        }
        else
        {
            if (DEBUGTEST) fprintf(stderr, "threaded_worker %p returns\n", varg);
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(& master->mutex);
#endif
            return;
        }
    }
}

slong _fmpz_mpoly_addmul_multi_threaded(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct * Blist,
    const slong * Blengths,
    const slong Bnumseq,
    const slong Btotallen,
    const flint_bitcnt_t bits,
    const slong N,
    const ulong * cmpmask,
    const fmpz_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_threads)
{
    struct _fmpz_mpoly_addmul_multi_worker * args;
    struct _fmpz_mpoly_addmul_multi_master * master;
    struct _fmpz_mpoly_addmul_multi_worker self;
    slong i;

    master = (struct _fmpz_mpoly_addmul_multi_master *) flint_malloc(sizeof(struct _fmpz_mpoly_addmul_multi_master));

    master->A = A;
    master->Blist = Blist;
    master->Blengths = Blengths;
    master->numterms = Bnumseq;
    master->Btotallen = Btotallen;
    master->bits = bits;
    master->N = N;
    master->cmpmask = cmpmask;
    master->ctx = ctx;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(& master->mutex, NULL);
#endif
    master->k = -WORD(1);
    master->control = (struct _fmpz_mpoly_addmul_multi_control *) flint_calloc(Bnumseq, sizeof(struct _fmpz_mpoly_addmul_multi_control));
    master->heap = NULL;
    master->merge_thread = NULL;

    args = (struct _fmpz_mpoly_addmul_multi_worker *) flint_malloc(num_threads
                                                  * sizeof(struct _fmpz_mpoly_addmul_multi_worker));
    for (i = 0; i < num_threads; i ++)
    {
        args[i].master = master;
        thread_pool_wake(global_thread_pool, handles[i], 0, _fmpz_mpoly_addmul_multi_threaded_worker, &args[i]);
    }

    self.master = master;
    _fmpz_mpoly_addmul_multi_threaded_worker(&self);

    for (i = 0; i < num_threads; i ++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    return master->k + 1;
}

void _fmpz_mpoly_addmul_multi_threaded_maxfields(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct ** Blist,
    const slong * Blengths,
    const slong Bnumseq,
    const slong Btotallen,
    const fmpz * maxfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N, Alen;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    int * freeBexp;
    ulong maxBlen = 0;
    slong maxlen = 0;
    int maxBindex;
    int aliasing_required = 0;
    slong i, j, k, m;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = Bnumseq;

    fmpz_mpoly_struct * B1;

    TMP_INIT;

    TMP_START;

    B1 = (fmpz_mpoly_struct *) TMP_ALLOC(Btotallen * sizeof(fmpz_mpoly_struct));

    Abits = _fmpz_vec_max_bits(maxfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);

    for (i = 0; i < Btotallen; i++)
    {
        Abits = FLINT_MAX(Abits, Blist[i]->bits);
        maxlen = maxlen + Blist[i]->length;

        if (Blist[i] == A)
            aliasing_required = 1;
    }

    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    freeBexp = (int *) TMP_ALLOC(Btotallen*sizeof(int));

    /* ensure input exponents are packed into same sized fields as output */
    k = 0;
    for (i = 0; i < Bnumseq; i++)
    {
        /* algorithm more efficient if largest poly first */
       maxBlen = 0;
       maxBindex = 0;
       for (j = 0; j < Blengths[i]; j++)
       {
           if (Blist[k + j]->length > maxBlen)
           {
               maxBlen = Blist[k + j]->length;
               maxBindex = k + j;
           }
       }

       for (j = 0; j < Blengths[i]; j++)
       {
           if (j == 0)
               m = maxBindex;
           else if (k + j <= maxBindex)
               m = k + j - 1;
           else
               m = k + j;

           B1[k + j].coeffs = Blist[m]->coeffs;
           B1[k + j].length = Blist[m]->length;
           B1[k + j].bits = Abits;

           if (Abits > Blist[m]->bits)
           {
               freeBexp[k + j] = 1;
               B1[k + j].exps = (ulong *) flint_malloc(N*Blist[m]->length*sizeof(ulong));
               mpoly_repack_monomials(B1[k + j].exps, Abits, Blist[m]->exps, Blist[m]->bits,
                                                               Blist[m]->length, ctx->minfo);
           }
           else
           {
               freeBexp[k + j] = 0;
               B1[k + j].exps = Blist[m]->exps;
           }
        }

        k += Blengths[i];
    }

    /* deal with aliasing and do multiplication */
    if (aliasing_required)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, maxlen, Abits, ctx);

        num_handles = flint_request_threads(&handles, thread_limit);
        if (DEBUGTEST) fprintf(stderr, "thread_limit %ld num_handles %ld\n", thread_limit, num_handles);
        Alen = _fmpz_mpoly_addmul_multi_threaded(T, B1, Blengths, Bnumseq, Btotallen, Abits, N, cmpmask, ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, maxlen, Abits, ctx);

        num_handles = flint_request_threads(&handles, thread_limit);
        if (DEBUGTEST) fprintf(stderr, "thread_limit %ld num_handles %ld\n", thread_limit, num_handles);
        Alen = _fmpz_mpoly_addmul_multi_threaded(A, B1, Blengths, Bnumseq, Btotallen, Abits, N, cmpmask, ctx, handles, num_handles);
        flint_give_back_threads(handles, num_handles);
    }

    for (i = 0; i < Btotallen; i++)
    {
        if (freeBexp[i])
            flint_free(B1[i].exps);
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}



/* B is a list of Blen polynomials to be multiplied together and the result is placed in A.
 *
 * Blist is an array of fmpz_mpoly_struct *, with Bnumseq subsequences of lengths given
 * by the lengths in the Blengths array.
 */

void fmpz_mpoly_addmul_multi_threaded(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct ** Blist,
    const slong * Blengths,
    const slong Bnumseq,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong Btotallen;
    fmpz * maxBfields;
    fmpz * maxtermfields;
    fmpz * maxfields;
    TMP_INIT;

    FLINT_ASSERT(Bnumseq > 0);

    if ((Bnumseq == 1) && (Blengths[0] == 1))
    {
        fmpz_mpoly_set(A, Blist[0], ctx);
	return;
    }

    TMP_START;

    flint_set_num_threads(2);

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxtermfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxtermfields + i);
        fmpz_init(maxfields + i);
    }

    /* Compute the maximum size of the exponent fields by computing the maximum size
     * of the fields for each term (by adding the max field sizes for each constituent polynomial),
     * then taking the max of all of them.
     */

    Btotallen = 0;
    k = 0;
    for (i = 0; i < Bnumseq; i++)
    {
        for (j = 0; j < Blengths[i]; j++)
        {
            mpoly_max_fields_fmpz(maxBfields, Blist[k]->exps, Blist[k]->length, Blist[k]->bits, ctx->minfo);
            _fmpz_vec_add(maxtermfields, maxtermfields, maxBfields, ctx->minfo->nfields);
            Btotallen ++;
            k ++;
        }
        _fmpz_vec_max(maxfields, maxfields, maxtermfields, ctx->minfo->nfields);
        _fmpz_vec_zero(maxtermfields, ctx->minfo->nfields);
    }

    _fmpz_mpoly_addmul_multi_threaded_maxfields(A, Blist, Blengths, Bnumseq, Btotallen, maxfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxtermfields + i);
        fmpz_clear(maxfields + i);
    }

    TMP_END;
}

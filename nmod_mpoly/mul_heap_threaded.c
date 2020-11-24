/*
    Copyright (C) 2017-2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <gmp.h>
#include <stdlib.h>
#include "thread_pool.h"
#include "nmod_mpoly.h"
#include "long_extras.h"


/*
    set A = the part of B*C with exps in [start, end)
    this functions reallocates A and returns the length of A
    version for N == 1
*/
void _nmod_mpoly_mul_heap_part1(
    nmod_mpoly_t A,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
    slong * start,
    slong * end,
    slong * hind,
    const nmod_mpoly_stripe_t S)
{
    const ulong cmpmask = S->cmpmask[0];
    slong i, j;
    ulong exp;
    mpoly_heap_t * x;
    slong next_loc;
    slong heap_len;
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = A->exps;
    mp_limb_t * Acoeff = A->coeffs;
    mp_limb_t acc0, acc1, acc2, pp0, pp1;

    FLINT_ASSERT(S->N == 1);

    /* tmp allocs from S->big_mem */
    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    heap = (mpoly_heap1_s *)(S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap1_s);
    chain = (mpoly_heap_t *)(S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    for (i = 0; i < Blen; i++)
    {
        hind[i] = 2*start[i] + 1;
    }
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;
      
        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, 1, Alen + 1);

        Aexp[Alen] = exp;

        acc0 = acc1 = acc2 = 0;
        do
        {
            x = _mpoly_heap_pop1(heap, &heap_len, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && heap[1].exp == exp);

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->mod);
        Alen += (Acoeff[Alen] != UWORD(0));

        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;
                hind[x->i] = 2*(x->j+1) + 0;
                _mpoly_heap_insert1(heap, Bexp[x->i] + Cexp[x->j], x,
                                                &next_loc, &heap_len, cmpmask);
            }
        }
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;
}


void _nmod_mpoly_mul_heap_part(
    nmod_mpoly_t A,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
    slong * start,
    slong * end,
    slong * hind,
    const nmod_mpoly_stripe_t S)
{
    flint_bitcnt_t bits = S->bits;
    slong N = S->N;
    const ulong * cmpmask = S->cmpmask;
    slong i, j;
    slong next_loc;
    slong heap_len;
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    mpoly_heap_t * x;
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    slong Alen;
    ulong * Aexp = A->exps;
    mp_limb_t * Acoeff = A->coeffs;
    ulong acc0, acc1, acc2, pp0, pp1;

    /* tmp allocs from S->big_mem */
    i = 0;
    store = store_base = (slong *) (S->big_mem + i);
    i += 2*Blen*sizeof(slong);
    exp_list = (ulong **) (S->big_mem + i);
    i += Blen*sizeof(ulong *);
    exps = (ulong *) (S->big_mem + i);
    i += Blen*N*sizeof(ulong);
    heap = (mpoly_heap_s *) (S->big_mem + i);
    i += (Blen + 1)*sizeof(mpoly_heap_s);
    chain = (mpoly_heap_t *) (S->big_mem + i);
    i += Blen*sizeof(mpoly_heap_t);
    FLINT_ASSERT(i <= S->big_mem_alloc);

    /* put all the starting nodes on the heap */
    heap_len = 1; /* heap zero index unused */
    next_loc = Blen + 4;   /* something bigger than heap can ever be */
    exp_next = 0;
    for (i = 0; i < Blen; i++)
        exp_list[i] = exps + N*i;
    for (i = 0; i < Blen; i++)
        hind[i] = 2*start[i] + 1;
    for (i = 0; i < Blen; i++)
    {
        if (  (start[i] < end[i])
           && (  (i == 0)
              || (start[i] < start[i - 1])
              )
           )
        {
            x = chain + i;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                       Cexp + N*x->j, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                          Cexp + N*x->j, N);

            exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
        }
    }

    Alen = 0;
    while (heap_len > 1)
    {
        exp = heap[1].exp;

        _nmod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, N, Alen + 1);

        mpoly_monomial_set(Aexp + N*Alen, exp, N);

        acc0 = acc1 = acc2 = 0;
        do
        {
            exp_list[--exp_next] = heap[1].exp;

            x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

            hind[x->i] |= WORD(1);
            *store++ = x->i;
            *store++ = x->j;
            umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
            add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);

            while ((x = x->next) != NULL)
            {
                hind[x->i] |= WORD(1);
                *store++ = x->i;
                *store++ = x->j;
                umul_ppmm(pp1, pp0, Bcoeff[x->i], Ccoeff[x->j]);
                add_sssaaaaaa(acc2, acc1, acc0, acc2, acc1, acc0, WORD(0), pp1, pp0);
            }
        } while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N));

        NMOD_RED3(Acoeff[Alen], acc2, acc1, acc0, S->mod);
        Alen += (Acoeff[Alen] != UWORD(0));

        /* for each node temporarily stored */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            /* should we go right? */
            if (  (i + 1 < Blen)
               && (j + 0 < end[i + 1])
               && (hind[i + 1] == 2*j + 1)
               )
            {
                x = chain + i + 1;
                x->i = i + 1;
                x->j = j;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }

            /* should we go up? */
            if (  (j + 1 < end[i + 0])
               && ((hind[i] & 1) == 1)
               && (  (i == 0)
                  || (hind[i - 1] >= 2*(j + 2) + 1)
                  )
               )
            {
                x = chain + i;
                x->i = i;
                x->j = j + 1;
                x->next = NULL;

                hind[x->i] = 2*(x->j + 1) + 0;

                if (bits <= FLINT_BITS)
                    mpoly_monomial_add(exp_list[exp_next], Bexp + N*x->i,
                                                           Cexp + N*x->j, N);
                else
                    mpoly_monomial_add_mp(exp_list[exp_next], Bexp + N*x->i,
                                                              Cexp + N*x->j, N);

                exp_next += _mpoly_heap_insert(heap, exp_list[exp_next], x,
                                             &next_loc, &heap_len, N, cmpmask);
            }
        }
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;
}


/*
    The workers calculate product terms from 4*n divisions, where n is the
    number of threads.
*/

typedef struct
{
    volatile int idx;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong nthreads;
    slong ndivs;
    const nmod_mpoly_ctx_struct * ctx;
    mp_limb_t * Acoeff;
    ulong * Aexp;
    const mp_limb_t * Bcoeff;
    const ulong * Bexp;
    slong Blen;
    const mp_limb_t * Ccoeff;
    const ulong * Cexp;
    slong Clen;
    slong N;
    flint_bitcnt_t bits;
    const ulong * cmpmask;
}
_base_struct;

typedef _base_struct _base_t[1];

typedef struct
{
    slong lower;
    slong upper;
    slong thread_idx;
    slong Aoffset;
    nmod_mpoly_t A;
}
_div_struct;

typedef struct
{
    nmod_mpoly_stripe_t S;
    slong idx;
    slong time;
    _base_struct * base;
    _div_struct * divs;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
    pthread_cond_t cond;
#endif
    slong * t1, * t2, * t3, * t4;
    ulong * exp;
}
_worker_arg_struct;


/*
    The workers simply take the next available division and calculate all
    product terms in this division.
*/

#define SWAP_PTRS(xx, yy) \
   do { \
      tt = xx; \
      xx = yy; \
      yy = tt; \
   } while (0)

static void _nmod_mpoly_mul_heap_threaded_worker(void * arg_ptr)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) arg_ptr;
    nmod_mpoly_stripe_struct * S = arg->S;
    _div_struct * divs = arg->divs;
    _base_struct * base = arg->base;
    slong Blen = base->Blen;
    slong N = base->N;
    slong i, j;
    ulong * exp;
    slong score;
    slong * start, * end, * t1, * t2, * t3, * t4, * tt;

    exp = (ulong *) flint_malloc(N*sizeof(ulong));
    t1 = (slong *) flint_malloc(Blen*sizeof(slong));
    t2 = (slong *) flint_malloc(Blen*sizeof(slong));
    t3 = (slong *) flint_malloc(Blen*sizeof(slong));
    t4 = (slong *) flint_malloc(Blen*sizeof(slong));

    S->N = N;
    S->bits = base->bits;
    S->cmpmask = base->cmpmask;
    S->ctx = base->ctx;
    S->mod = base->ctx->mod;

    S->big_mem_alloc = 0;
    if (N == 1)
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap1_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
    }
    else
    {
        S->big_mem_alloc += 2*Blen*sizeof(slong);
        S->big_mem_alloc += (Blen + 1)*sizeof(mpoly_heap_s);
        S->big_mem_alloc += Blen*sizeof(mpoly_heap_t);
        S->big_mem_alloc += Blen*S->N*sizeof(ulong);
        S->big_mem_alloc += Blen*sizeof(ulong *);
    }
    S->big_mem = (char *) flint_malloc(S->big_mem_alloc);

    /* get index to start working on */
    if (arg->idx + 1 < base->nthreads)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }
    else
    {
        i = base->ndivs - 1;
    }

    while (i >= 0)
    {
        FLINT_ASSERT(divs[i].thread_idx == -WORD(1));
        divs[i].thread_idx = arg->idx;

        /* calculate start */
        if (i + 1 < base-> ndivs)
        {
            mpoly_search_monomials(
                &start, exp, &score, t1, t2, t3,
                            divs[i].lower, divs[i].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
                                          base->N, base->cmpmask);
            if (start == t2)
            {
                SWAP_PTRS(t1, t2);
            }
            else if (start == t3)
            {
                SWAP_PTRS(t1, t3);
            }
        }
        else
        {
            start = t1;
            for (j = 0; j < base->Blen; j++)
                start[j] = 0;
        }

        /* calculate end */
        if (i > 0)
        {
            mpoly_search_monomials(
                &end, exp, &score, t2, t3, t4,
                            divs[i - 1].lower, divs[i - 1].lower,
                            base->Bexp, base->Blen, base->Cexp, base->Clen,
                                          base->N, base->cmpmask);


            if (end == t3)
            {
                SWAP_PTRS(t2, t3);
            }
            else if (end == t4)
            {
                SWAP_PTRS(t2, t4);
            }
        }
        else
        {
            end = t2;
            for (j = 0; j < base->Blen; j++)
                end[j] = base->Clen;
        }
        /* t3 and t4 are free for workspace at this point */

        /* calculate products in [start, end) */
        _nmod_mpoly_fit_length(&divs[i].A->coeffs, &divs[i].A->coeffs_alloc,
                             &divs[i].A->exps, &divs[i].A->exps_alloc, N, 256);
        if (N == 1)
        {
            _nmod_mpoly_mul_heap_part1(divs[i].A,
                                      base->Bcoeff,  base->Bexp,  base->Blen,
                                      base->Ccoeff,  base->Cexp,  base->Clen,
                                                            start, end, t3, S);
        }
        else
        {
            _nmod_mpoly_mul_heap_part(divs[i].A,
                                      base->Bcoeff,  base->Bexp,  base->Blen,
                                      base->Ccoeff,  base->Cexp,  base->Clen,
                                                            start, end, t3, S);
        }

        /* get next index to work on */
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&base->mutex);
#endif
	i = base->idx - 1;
        base->idx = i;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&base->mutex);
#endif
    }

    /* clean up */
    flint_free(S->big_mem);
    flint_free(t4);
    flint_free(t3);
    flint_free(t2);
    flint_free(t1);
    flint_free(exp);
}


static void _join_worker(void * varg)
{
    _worker_arg_struct * arg = (_worker_arg_struct *) varg;
    _div_struct * divs = arg->divs;
    _base_struct * base = arg->base;
    slong N = base->N;
    slong i;

    for (i = base->ndivs - 2; i >= 0; i--)
    {
        FLINT_ASSERT(divs[i].thread_idx != -WORD(1));

        if (divs[i].thread_idx != arg->idx)
            continue;

        FLINT_ASSERT(divs[i].A->coeffs != NULL);
        FLINT_ASSERT(divs[i].A->exps != NULL);

        memcpy(base->Acoeff + divs[i].Aoffset, divs[i].A->coeffs,
                                          divs[i].A->length*sizeof(mp_limb_t));

        memcpy(base->Aexp + N*divs[i].Aoffset, divs[i].A->exps,
                                            N*divs[i].A->length*sizeof(ulong));

        flint_free(divs[i].A->coeffs);
        flint_free(divs[i].A->exps);
    }
}

void _nmod_mpoly_mul_heap_threaded(
    nmod_mpoly_t A,
    const mp_limb_t * Bcoeff, const ulong * Bexp, slong Blen,
    const mp_limb_t * Ccoeff, const ulong * Cexp, slong Clen,
    flint_bitcnt_t bits,
    slong N,
    const ulong * cmpmask,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong i;
    _base_t base;
    _div_struct * divs;
    _worker_arg_struct * args;
    slong BClen, Alen;

    /* bail if product of lengths overflows a word */
    if (z_mul_checked(&BClen, Blen, Clen))
    {
        _nmod_mpoly_mul_johnson(A, Bcoeff, Bexp, Blen,
                               Ccoeff, Cexp, Clen, bits, N, cmpmask, ctx->mod);
        return;
    }

    base->nthreads = num_handles + 1;
    base->ndivs    = base->nthreads*4;  /* number of divisons */
    base->Bcoeff = Bcoeff;
    base->Bexp = Bexp;
    base->Blen = Blen;
    base->Ccoeff = Ccoeff;
    base->Cexp = Cexp;
    base->Clen = Clen;
    base->bits = bits;
    base->N = N;
    base->cmpmask = cmpmask;
    base->idx = base->ndivs - 1;    /* decremented by worker threads */
    base->ctx = ctx;

    divs = (_div_struct *) flint_malloc(base->ndivs*sizeof(_div_struct));
    args = (_worker_arg_struct *) flint_malloc(base->nthreads
                                                  *sizeof(_worker_arg_struct));

    /* allocate space and set the boundary for each division */
    FLINT_ASSERT(BClen/Blen == Clen);
    for (i = base->ndivs - 1; i >= 0; i--)
    {
        double d = (double)(i + 1) / (double)(base->ndivs);

        /* divisions decrease in size so that no worker finishes too early */
        divs[i].lower = (d * d) * BClen;
        divs[i].lower = FLINT_MIN(divs[i].lower, BClen);
        divs[i].lower = FLINT_MAX(divs[i].lower, WORD(0));
        divs[i].upper = divs[i].lower;
        divs[i].Aoffset = -WORD(1);
        divs[i].thread_idx = -WORD(1);

        if (i == base->ndivs - 1)
        {
            /* highest division writes to original poly */
            *divs[i].A = *A;
        }
        else
        {
            /* lower divisions write to a new worker poly */
            divs[i].A->exps = NULL;
            divs[i].A->coeffs = NULL;
            divs[i].A->bits = A->bits;
            divs[i].A->coeffs_alloc = 0;
            divs[i].A->exps_alloc = 0;
        }
        divs[i].A->length = 0;
    }

    /* compute each chunk in parallel */
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&base->mutex, NULL);
#endif
    for (i = 0; i < num_handles; i++)
    {
        args[i].idx = i;
        args[i].base = base;
        args[i].divs = divs;
        thread_pool_wake(global_thread_pool, handles[i], 0,
                               _nmod_mpoly_mul_heap_threaded_worker, &args[i]);
    }
    i = num_handles;
    args[i].idx = i;
    args[i].base = base;
    args[i].divs = divs;
    _nmod_mpoly_mul_heap_threaded_worker(&args[i]);
    for (i = 0; i < num_handles; i++)
    {
        thread_pool_wait(global_thread_pool, handles[i]);
    }

    /* calculate and allocate space for final answer */
    i = base->ndivs - 1;
    *A = *divs[i].A;
    Alen = A->length;
    for (i = base->ndivs - 2; i >= 0; i--)
    {
        divs[i].Aoffset = Alen;
        Alen += divs[i].A->length;
    }
    nmod_mpoly_fit_length(A, Alen, ctx);

    base->Acoeff = A->coeffs;
    base->Aexp = A->exps;

    /* join answers */
    for (i = 0; i < num_handles; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, _join_worker, &args[i]);
    _join_worker(&args[num_handles]);
    for (i = 0; i < num_handles; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    A->length = Alen;

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&base->mutex);
#endif

    flint_free(args);
    flint_free(divs);
}


/* maxBfields gets clobbered */
void _nmod_mpoly_mul_heap_threaded_pool_maxfields(
    nmod_mpoly_t A,
    const nmod_mpoly_t B, fmpz * maxBfields,
    const nmod_mpoly_t C, fmpz * maxCfields,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    slong N;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    ulong * Bexp, * Cexp;
    int freeBexp, freeCexp;
    nmod_mpoly_struct * a, T[1];

    TMP_INIT;

    TMP_START;

    _fmpz_vec_add(maxBfields, maxBfields, maxCfields, ctx->minfo->nfields);

    Abits = 1 + _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(Abits, B->bits);
    Abits = FLINT_MAX(Abits, C->bits);
    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    /* ensure input exponents are packed into same sized fields as output */
    freeBexp = 0;
    Bexp = B->exps;
    if (Abits > B->bits)
    {
        freeBexp = 1;
        Bexp = (ulong *) flint_malloc(N*B->length*sizeof(ulong));
        mpoly_repack_monomials(Bexp, Abits, B->exps, B->bits, B->length, ctx->minfo);
    }

    freeCexp = 0;
    Cexp = C->exps;
    if (Abits > C->bits)
    {
        freeCexp = 1;
        Cexp = (ulong *) flint_malloc(N*C->length*sizeof(ulong));
        mpoly_repack_monomials(Cexp, Abits, C->exps, C->bits, C->length, ctx->minfo);
    }

    if (A == B || A == C)
    {
        nmod_mpoly_init(T, ctx);
        a = T;
    }
    else
    {
        a = A;
    }

    nmod_mpoly_fit_length_reset_bits(a, B->length + C->length, Abits, ctx);

    if (B->length > C->length)
    {
        _nmod_mpoly_mul_heap_threaded(a, C->coeffs, Cexp, C->length,
                                         B->coeffs, Bexp, B->length,
                             Abits, N, cmpmask, ctx, handles, num_handles);
    }
    else
    {
        _nmod_mpoly_mul_heap_threaded(a, B->coeffs, Bexp, B->length,
                                         C->coeffs, Cexp, C->length,
                             Abits, N, cmpmask, ctx, handles, num_handles);
    }

    if (A == B || A == C)
    {
        nmod_mpoly_swap(A, T, ctx);
        nmod_mpoly_clear(T, ctx);
    }

    if (freeBexp)
        flint_free(Bexp);

    if (freeCexp)
        flint_free(Cexp);

    TMP_END;
}


void nmod_mpoly_mul_heap_threaded(
    nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_t C,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxBfields, * maxCfields;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(B->length, C->length)/16;
    TMP_INIT;

    if (B->length == 0 || C->length == 0)
    {
        nmod_mpoly_zero(A, ctx);
        return;
    }

    TMP_START;

    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);

    num_handles = flint_request_threads(&handles, thread_limit);

    _nmod_mpoly_mul_heap_threaded_pool_maxfields(A, B, maxBfields, C, maxCfields,
                                                    ctx, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }

    TMP_END;
}

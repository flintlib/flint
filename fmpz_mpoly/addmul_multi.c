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
*/
slong _fmpz_mpoly_addmul_multi(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct * Blist,
    const slong * Blengths,
    const slong Bnumseq,
    const slong Btotallen,
    const flint_bitcnt_t bits,
    const slong N,
    const ulong * cmpmask,
    const fmpz_mpoly_ctx_t ctx)
{
   slong i, j, k, l;
   slong next_loc;
   slong Q_len = 0, heap_len = 1; /* heap starts empty, and its zero index is unused, so heap_len = 1 */
   slong heap_size;
   slong heap_block_size;
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** chain_list;
   ulong * Q;
   slong Q_size;
   mpoly_heap_t * x;
   ulong multiindex;
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   slong chain_next;
   slong chain_size;
   hind_t * hind;
   int first;
   slong hind_totallen = 0;
   slong * hind_len;
   slong * hind_start;
   ulong offset, offset2, offset3;
   ulong candidate;
   ulong partial_multiindex;
   fmpz_t tmp_coeff;
   const fmpz_mpoly_struct ** Bstart;
   slong term;
   TMP_INIT;

   TMP_START;

   fmpz_init(tmp_coeff);

   Bstart = (const fmpz_mpoly_struct **) TMP_ALLOC(Bnumseq*sizeof(fmpz_mpoly_struct *));

   /* First polynomial should be the largest; it's left out of this calculation */
   hind_len = TMP_ALLOC(Bnumseq*sizeof(slong));
   hind_start = TMP_ALLOC(Bnumseq*sizeof(slong));

   for (i = 0; i < Bnumseq; i++)
   {
      if (i == 0)
          Bstart[i] = Blist;
      else
          Bstart[i] = Bstart[i-1] + Blengths[i-1];

      hind_start[i] = hind_totallen;
      hind_len[i] = 1;

      /* XXX - FLINT_ASSERTs are only tested if configured with --enable-assert */
      FLINT_ASSERT(Bstart[i][0].length < (UWORD(1) << (8*sizeof(hind_t) - 1)));

      for (j = 1; j < Blengths[i]; j++)
      {
          FLINT_ASSERT(Bstart[i][j].length <= Bstart[i][0].length);
          hind_len[i] *= Bstart[i][j].length;
      }
      hind_totallen += hind_len[i];
   }

   next_loc = hind_totallen + 4;   /* something bigger than heap can ever be */

   heap_block_size = 16384/sizeof(ulong)/N;
   for (heap_size = 0; heap_size < Bnumseq + 1; heap_size += heap_block_size);
   for (chain_size = 0; chain_size < Bnumseq; chain_size += chain_block_size);

   heap = (mpoly_heap_s *) flint_malloc(heap_size*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(chain_size*sizeof(mpoly_heap_t));
   chain_list = (mpoly_heap_t **) flint_malloc(chain_size*sizeof(mpoly_heap_t **));
   for (i = 0; i < chain_size; i++)
      chain_list[i] = chain + i;

   /* space for temporary storage of pointers to heap nodes */
   Q_size = Bnumseq;
   Q = (ulong *) flint_malloc(Q_size*sizeof(ulong));

   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(heap_size*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) flint_malloc(heap_size*sizeof(ulong *));
   for (i = 0; i < heap_size; i++)
      exp_list[i] = exps + i*N;

   /* space for heap indices */
   hind = (hind_t *) TMP_ALLOC(hind_totallen*sizeof(hind_t));
   for (i = 0; i < hind_totallen; i++)
       hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   exp_next = 0;
   chain_next = 0;

   for (i = 0; i < Bnumseq; i++)
   {
      /* for each term, put (0...0, exp2[0] + ... + expn[0]) on heap */
      x = chain_list[chain_next ++];
      x->i = hind_start[i];
      x->next = NULL;

      mpoly_monomial_zero(exp_list[exp_next], N);
      for (j = 0; j < Blengths[i]; j++)
          if (bits <= FLINT_BITS)
              mpoly_monomial_add(exp_list[exp_next], exp_list[exp_next], Bstart[i][j].exps, N);
          else
              mpoly_monomial_add_mp(exp_list[exp_next], exp_list[exp_next], Bstart[i][j].exps, N);

      hind[hind_start[i] + 0] = 2*1 + 0;

      if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                              &next_loc, &heap_len, N, cmpmask))
         exp_next--;
   }

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get pointer to exponent field of heap top */
      exp = heap[1].exp;

      /* realloc output poly ready for next product term */
      k++;
      fmpz_mpoly_fit_length(A, k + 1, ctx);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* pop chain from heap and set exponent field to be reused */
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N, cmpmask);

         if (first)
         {
            fmpz_zero(A->coeffs + k);
            
            /* set output monomial */
            mpoly_monomial_set(A->exps + k*N, exp, N);

            first = 0; 
         }

         /* for each node in this chain */
         do
         {
            multiindex = x->i;

            for (term=0; term < Bnumseq-1; term++)
            {
                if (multiindex < hind_start[term+1])
                    break;
            }

            if (EXPTEST) {
                fprintf(stderr, "Processing %ld from term %ld with exp %lx\n", multiindex, term, *exp);
            }
            /* addmul product of input poly coeffs */

            fmpz_set(tmp_coeff, Bstart[term][0].coeffs + (hind[multiindex] >> 1) - 1);
            partial_multiindex = multiindex - hind_start[term];
            for (i=1; i<Blengths[term]; i++)
            {
                fmpz_mul(tmp_coeff, tmp_coeff, Bstart[term][i].coeffs + (partial_multiindex % Bstart[term][i].length));
                partial_multiindex /= Bstart[term][i].length;
            }
            fmpz_add(A->coeffs + k, A->coeffs + k, tmp_coeff);

            if (EXPTEST) {
                fprintf(stderr, "adding tmp_coeff = ");
                fmpz_fprint(stderr, tmp_coeff);
                fprintf(stderr, " to get A->coeffs[%ld]=", k);
                fmpz_fprint(stderr, A->coeffs + k);
                fprintf(stderr, "\n");
            }

            /* take node out of heap and put into store */
            if (Q_len == Q_size)
            {
                Q_size += Bnumseq;
                Q = (ulong *) flint_realloc(Q, Q_size*sizeof(ulong));
            }
            hind[multiindex] |= WORD(1);
            Q[Q_len++] = multiindex;

            chain_list[-- chain_next] = x;
         }
         while ((x = x->next) != NULL);
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         multiindex = Q[--Q_len];

         for (term=0; term < Bnumseq-1; term++)
         {
             if (multiindex < hind_start[term+1])
                 break;
         }

         offset = 0;
         for (i=0; i<Blengths[term]; i++)
         {
             candidate = multiindex - hind_start[term] + offset;
             if (EXPTEST) fprintf(stderr, "Considering term %ld candidate %ld\n", term, candidate);
             if ((candidate < hind_len[term]) && (hind[hind_start[term] + candidate] < 2*Bstart[term][0].length) && (hind[hind_start[term] + candidate] & 1))
             {
                 offset2 = 1;
                 for (j=1; j<Blengths[term]; j++)
                 {
                     offset3 = offset2 * Bstart[term][j].length;
                     if (((candidate % offset3) / offset2 != 0) && (hind[hind_start[term] + candidate - offset2] < hind[hind_start[term] + candidate] + 2))
                         break;
                     offset2 *= Bstart[term][j].length;
                 }
                 if (j == Blengths[term])
                 {
                     if (heap_len == heap_size)
                     {
                         heap_size += heap_block_size;
                         heap = flint_realloc(heap, heap_size*sizeof(mpoly_heap_s));
                         exps = TMP_ALLOC(heap_block_size*N*sizeof(ulong));
                         exp_list = flint_realloc(exp_list, heap_size*sizeof(ulong *));
                         for (l = 0; l < heap_block_size; l++)
                             exp_list[heap_len + l] = exps + l*N;
                     }

                     if (chain_next == chain_size)
                     {
                         chain_size += chain_block_size;
                         chain = (mpoly_heap_t *) TMP_ALLOC(chain_block_size*sizeof(mpoly_heap_t));
                         chain_list = (mpoly_heap_t **) flint_realloc(chain_list, chain_size*sizeof(mpoly_heap_t *));
                         for (l = 0; l < chain_block_size; l++)
                             chain_list[chain_next + l] = chain + l;
                     }

                     x = chain_list[chain_next ++];
                     x->i = hind_start[term] + candidate;
                     x->next = NULL;

                     mpoly_monomial_set(exp_list[exp_next], Bstart[term][0].exps + N*(hind[hind_start[term] + candidate] >> 1), N);
                     hind[hind_start[term] + candidate] ++;

                     partial_multiindex = candidate;
                     for (l = 1; l < Blengths[term]; l++)
                     {
                         if (bits <= FLINT_BITS)
                             mpoly_monomial_add(exp_list[exp_next], exp_list[exp_next],
                                                Bstart[term][l].exps + N*(partial_multiindex % Bstart[term][l].length), N);
                         else
                             mpoly_monomial_add_mp(exp_list[exp_next], exp_list[exp_next],
                                                   Bstart[term][l].exps + N*(partial_multiindex % Bstart[term][l].length), N);
                         partial_multiindex /= Bstart[term][l].length;
                     }

                     if (EXPTEST) {
	                 fprintf(stderr, "Adding %ld (term %ld) because of %ld with hind[%ld] " HIND_FORMAT "\n",
                                 candidate, term, multiindex, candidate, hind[hind_start[term] + candidate]);
                     }
                     if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                             &next_loc, &heap_len, N, cmpmask))
                         exp_next--;
                 }
             }

             if (offset == 0)
                 offset = 1;
             else
                 offset *= Bstart[term][i].length;
         }
      }

      if (fmpz_is_zero(A->coeffs + k))
         k--;
   }

   k++;

   fmpz_clear(tmp_coeff);

   flint_free(exp_list);
   flint_free(Q);
   flint_free(chain_list);
   flint_free(heap);

   TMP_END;

   return k;
}

void _fmpz_mpoly_addmul_multi_maxfields(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct ** Blist,
    const slong * Blengths,
    const slong Bnumseq,
    slong Btotallen,
    const fmpz * maxfields,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N, Alen;
    flint_bitcnt_t Abits;
    ulong * cmpmask;
    int * freeBexp;
    ulong numterms;
    ulong maxBlen = 0;
    slong maxlen = 0;
    slong * B1lengths;
    int maxBindex;
    int aliasing_required = 0;
    slong i, j, k, k1, m;

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

    B1lengths = (slong *) TMP_ALLOC(Bnumseq * sizeof(slong));

    /* ensure input exponents are packed into same sized fields as output */
    k = 0;
    k1 = 0;
    numterms = 0;
    for (i = 0; i < Bnumseq; i++)
    {
        /* algorithm more efficient if largest poly first */
       maxBlen = 0;
       maxBindex = 0;
       for (j = 0; j < Blengths[i]; j++)
       {
           if (Blist[k + j]->length == 0)
               break;
           if (Blist[k + j]->length > maxBlen)
           {
               maxBlen = Blist[k + j]->length;
               maxBindex = k + j;
           }
       }

       /* FLINT's convention is that zero-length polynomials are treated as zero polynomials,
        * so drop any term that contained a zero-length polynomial as a factor.
        */

       if (j < Blengths[i])
       {
           k += Blengths[i];
           Btotallen -= Blengths[i];
           continue;
       }

       for (j = 0; j < Blengths[i]; j++)
       {
           if (j == 0)
               m = maxBindex;
           else if (k + j <= maxBindex)
               m = k + j - 1;
           else
               m = k + j;

           B1[k1 + j].coeffs = Blist[m]->coeffs;
           B1[k1 + j].length = Blist[m]->length;
           B1[k1 + j].bits = Abits;

           if (Abits > Blist[m]->bits)
           {
               freeBexp[k1 + j] = 1;
               B1[k1 + j].exps = (ulong *) flint_malloc(N*Blist[m]->length*sizeof(ulong));
               mpoly_repack_monomials(B1[k1 + j].exps, Abits, Blist[m]->exps, Blist[m]->bits,
                                                               Blist[m]->length, ctx->minfo);
           }
           else
           {
               freeBexp[k + j] = 0;
               B1[k1 + j].exps = Blist[m]->exps;
           }
        }

        k += Blengths[i];
        k1 += Blengths[i];
        B1lengths[numterms ++] = Blengths[i];
    }

    if (Btotallen == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    /* deal with aliasing and do multiplication */
    if (aliasing_required)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, maxlen, Abits, ctx);

        Alen = _fmpz_mpoly_addmul_multi(T, B1, B1lengths, numterms, Btotallen, Abits, N, cmpmask, ctx);

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, maxlen, Abits, ctx);

        Alen = _fmpz_mpoly_addmul_multi(A, B1, B1lengths, numterms, Btotallen, Abits, N, cmpmask, ctx);
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

void fmpz_mpoly_addmul_multi(
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

    _fmpz_mpoly_addmul_multi_maxfields(A, Blist, Blengths, Bnumseq, Btotallen, maxfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxtermfields + i);
        fmpz_clear(maxfields + i);
    }

    TMP_END;
}

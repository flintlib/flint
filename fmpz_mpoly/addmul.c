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

/*
   Set A to B1*B2*---*Bn using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors take N words.
*/
slong _fmpz_mpoly_addmul(
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
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   ulong * Q;
   mpoly_heap_t * x;
   ulong multiindex;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   slong * hind;
   int first;
   slong hind_len;
   ulong offset, offset2, offset3;
   ulong candidate;
   ulong partial_multiindex;
   fmpz_t tmp_coeff;
   TMP_INIT;

   TMP_START;

   /* for right now, we only can do a single multiplication */
   FLINT_ASSERT(Blengths[0] == Btotallen);

   fmpz_init(tmp_coeff);

   /* First polynomial should be the largest; it's left out of this calculation */
   hind_len = 1;
   for (i = 1; i < Btotallen; i++)
       hind_len *= Blist[i].length;

   next_loc = hind_len + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap_s *) TMP_ALLOC((hind_len + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(hind_len*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (ulong *) TMP_ALLOC(hind_len*sizeof(ulong));
   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(hind_len*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(hind_len*sizeof(ulong *));
   for (i = 0; i < hind_len; i++)
      exp_list[i] = exps + i*N;

   /* space for heap indices */
   hind = (slong *) TMP_ALLOC(hind_len*sizeof(slong));
   for (i = 0; i < hind_len; i++)
       hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   exp_next = 0;

   /* put (0...0, exp2[0] + ... + expn[0]) on heap */
   x = chain + 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

   mpoly_monomial_zero(heap[1].exp, N);
   for (i = 0; i < Btotallen; i++)
       if (bits <= FLINT_BITS)
           mpoly_monomial_add(heap[1].exp, heap[1].exp, Blist[i].exps, N);
       else
           mpoly_monomial_add_mp(heap[1].exp, heap[1].exp, Blist[i].exps, N);

   hind[0] = 2*1 + 0;

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

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

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
            multiindex = x - chain;

            if (EXPTEST) {
                fprintf(stderr, "Processing %ld with exp %lx\n", multiindex, *exp);
            }
            /* addmul product of input poly coeffs */

            fmpz_set(tmp_coeff, Blist[0].coeffs + (hind[multiindex] >> 1) - 1);
            partial_multiindex = multiindex;
            for (i=1; i<Btotallen-1; i++)
            {
                fmpz_mul(tmp_coeff, tmp_coeff, Blist[i].coeffs + (partial_multiindex % Blist[i].length));
                partial_multiindex /= Blist[i].length;
            }
            fmpz_mul(tmp_coeff, tmp_coeff, Blist[Btotallen-1].coeffs + partial_multiindex);
            fmpz_add(A->coeffs + k, A->coeffs + k, tmp_coeff);

            if (EXPTEST) {
                fprintf(stderr, "adding tmp_coeff = ");
                fmpz_fprint(stderr, tmp_coeff);
                fprintf(stderr, " to get A->coeffs[%ld]=", k);
                fmpz_fprint(stderr, A->coeffs + k);
                fprintf(stderr, "\n");
            }

            /* take node out of heap and put into store */
            hind[multiindex] |= WORD(1);
            Q[Q_len++] = multiindex;
         }
         while ((x = x->next) != NULL);
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         multiindex = Q[--Q_len];

         offset = 0;
         for (i=0; i<Btotallen; i++)
         {
             candidate = multiindex + offset;
             if (EXPTEST) fprintf(stderr, "Considering candidate %ld\n", candidate);
             if ((candidate < hind_len) && (hind[candidate] < 2*Blist[0].length) && (hind[candidate] & 1))
             {
                 offset2 = 1;
                 for (j=1; j<Btotallen; j++)
                 {
                     offset3 = offset2 * Blist[j].length;
                     if (((candidate % offset3) / offset2 != 0) && (hind[candidate - offset2] < hind[candidate] + 2))
                         break;
                     offset2 *= Blist[j].length;
                 }
                 if (j == Btotallen)
                 {
                     x = chain + candidate;
                     x->next = NULL;

                     mpoly_monomial_set(exp_list[exp_next], Blist[0].exps + N*(hind[candidate] >> 1), N);
                     hind[candidate] ++;

                     partial_multiindex = candidate;
                     for (l = 1; l < Btotallen; l++)
                     {
                         if (bits <= FLINT_BITS)
                             mpoly_monomial_add(exp_list[exp_next], exp_list[exp_next],
                                                Blist[l].exps + N*(partial_multiindex % Blist[l].length), N);
                         else
                             mpoly_monomial_add_mp(exp_list[exp_next], exp_list[exp_next],
                                                   Blist[l].exps + N*(partial_multiindex % Blist[l].length), N);
                         partial_multiindex /= Blist[l].length;
                     }

                     if (EXPTEST) {
	                 fprintf(stderr, "Adding %ld because of %ld with hind[%ld] %ld\n", candidate, multiindex, candidate, hind[candidate]);
                     }
                     if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                             &next_loc, &heap_len, N, cmpmask))
                         exp_next--;
                 }
             }

             if (offset == 0)
                 offset = 1;
             else
                 offset *= Blist[i].length;
         }
      }

      if (fmpz_is_zero(A->coeffs + k))
         k--;
   }

   k++;

   fmpz_clear(tmp_coeff);

   TMP_END;

   return k;
}

void _fmpz_mpoly_addmul_maxfields(
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
           if (k + j == maxBindex)
               m = k;
           else if (k + j < maxBindex)
               m = k + j + 1;
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

        Alen = _fmpz_mpoly_addmul(T, B1, Blengths, Bnumseq, Btotallen, Abits, N, cmpmask, ctx);

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, maxlen, Abits, ctx);

        Alen = _fmpz_mpoly_addmul(A, B1, Blengths, Bnumseq, Btotallen, Abits, N, cmpmask, ctx);
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

void fmpz_mpoly_addmul(
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

    _fmpz_mpoly_addmul_maxfields(A, Blist, Blengths, Bnumseq, Btotallen, maxfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxtermfields + i);
        fmpz_clear(maxfields + i);
    }

    TMP_END;
}

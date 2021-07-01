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
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors all fit in a
   single word. Assumes input polys are nonzero.
*/
slong _fmpz_mpoly_addmul1(fmpz ** poly1, ulong ** exp1, slong * alloc,
              const fmpz * poly2, const ulong * exp2, slong len2,
              const fmpz * poly3, const ulong * exp3, slong len3, ulong maskhi)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   slong * hind;
   ulong exp, cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, small;
   TMP_INIT;

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));

    /* space for heap indices */
    hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
    for (i = 0; i < len2; i++)
        hind[i] = 1;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0] + exp3[0], x);
   hind[0] = 2*1 + 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get exponent field of heap top */
      exp = heap[1].exp;
      
      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && heap[1].exp == exp)
      {
         /* pop chain from heap */
         x = _mpoly_heap_pop1(heap, &heap_len, maskhi);
         
         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               e1[k] = exp;
               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }

            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               e1[k] = exp;
               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }
      
      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         j = Q[--Q_len];
         i = Q[--Q_len];

         /* should we go right? */
         if (  (i + 1 < len2)
            && (hind[i + 1] == 2*j + 1)
            )
         {
            x = chain + i + 1;
            x->i = i + 1;
            x->j = j;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }

         /* should we go up? */
         if (  (j + 1 < len3)
            && ((hind[i] & 1) == 1)
            && (  (i == 0)
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
               )
            )
         {
            x = chain + i;
            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x,
                                                 &next_loc, &heap_len, maskhi);
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors take N words.
*/
slong _fmpz_mpoly_addmul(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct * B,
    const slong Blen,
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
   int first, small;
   slong hind_len;
   ulong offset, offset2;
   ulong candidate;
   ulong partial_multiindex;
   fmpz_t tmp_coeff;
   TMP_INIT;

   /* if exponent vectors fit in single word, call special version */
#if 0
   if (N == 1)
      return _fmpz_mpoly_addmul1(poly1, exp1, alloc,
                             poly2, exp2, len2, poly3, exp3, len3, cmpmask[0]);
#endif

   TMP_START;

   fmpz_init(tmp_coeff);

   /* whether input coeffs are small, thus output coeffs fit in three words */
#if 0
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);
#else
   small = 0;
#endif

   /* First polynomial should be the largest; it's left out of this calculation */
   hind_len = 1;
   for (i = 1; i < Blen; i++)
       hind_len *= B[i].length;

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
   for (i = 0; i < Blen; i++)
       if (bits <= FLINT_BITS)
           mpoly_monomial_add(heap[1].exp, heap[1].exp, B[i].exps, N);
       else
           mpoly_monomial_add_mp(heap[1].exp, heap[1].exp, B[i].exps, N);

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

         /* if output coeffs will fit in three words */
         if (small)
         {
#if 0
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               mpoly_monomial_set(A->exps + k*N, exp, N);

               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
#endif
         } else /* output coeffs require multiprecision */
         {
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

               fmpz_set(tmp_coeff, B[0].coeffs + (hind[multiindex] >> 1) - 1);
               partial_multiindex = multiindex;
               for (i=1; i<Blen-1; i++)
               {
                   fmpz_mul(tmp_coeff, tmp_coeff, B[i].coeffs + (partial_multiindex % B[i].length));
                   partial_multiindex /= B[i].length;
               }
               fmpz_mul(tmp_coeff, tmp_coeff, B[Blen-1].coeffs + partial_multiindex);
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
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         multiindex = Q[--Q_len];

         offset = 0;
         for (i=0; i<Blen; i++)
         {
             candidate = multiindex + offset;
             if ((candidate < hind_len) && (hind[candidate] < 2*B[0].length) && (hind[candidate] & 1))
             {
                 offset2 = 1;
                 for (j=1; j<Blen; j++)
                 {
                     if (((candidate - offset2) >= 0) && (hind[candidate - offset2] < hind[candidate] + 2))
                         break;
                     offset2 *= B[i].length;
                 }
                 if (j == Blen)
                 {
                     x = chain + candidate;
                     x->next = NULL;

                     mpoly_monomial_set(exp_list[exp_next], B[0].exps + (hind[candidate] >> 1), N);
                     hind[candidate] ++;

                     partial_multiindex = candidate;
                     for (l = 1; l < Blen; l++)
                     {
                         if (bits <= FLINT_BITS)
                             mpoly_monomial_add(exp_list[exp_next], exp_list[exp_next],
                                                B[l].exps + (partial_multiindex % B[l].length), N);
                         else
                             mpoly_monomial_add_mp(exp_list[exp_next], exp_list[exp_next],
                                                   B[l].exps + (partial_multiindex % B[l].length), N);
                         partial_multiindex /= B[l].length;
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
                 offset *= B[i].length;
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(A->coeffs + k, c[2], c[1], c[0]);

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
    const fmpz_mpoly_struct * B,
    const slong Blen,
    const fmpz * maxBfields,
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
    slong i;
    slong j;

    fmpz_mpoly_struct * B1;

    TMP_INIT;

    TMP_START;

    B1 = (fmpz_mpoly_struct *) TMP_ALLOC(Blen * sizeof(fmpz_mpoly_struct));

    Abits = _fmpz_vec_max_bits(maxBfields, ctx->minfo->nfields);
    Abits = FLINT_MAX(MPOLY_MIN_BITS, Abits + 1);

    for (i = 0; i < Blen; i++)
    {
        Abits = FLINT_MAX(Abits, B[i].bits);
        maxlen = maxlen + B[i].length;
        if (maxBlen < B[i].length)
        {
            maxBlen = B[i].length;
            maxBindex = i;
        }
        if ((B+i) == A)
            aliasing_required = 1;
    }

    Abits = mpoly_fix_bits(Abits, ctx->minfo);

    N = mpoly_words_per_exp(Abits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, Abits, ctx->minfo);

    freeBexp = (int *) TMP_ALLOC(Blen*sizeof(int));

    /* ensure input exponents are packed into same sized fields as output */
    for (i = 0; i < Blen; i++)
    {
        /* algorithm more efficient if largest poly first */
       if (i == maxBindex)
           j = 0;
       else if (i < maxBindex)
           j = i + 1;
       else
           j = i;

       B1[j].coeffs = B[i].coeffs;
       B1[j].length = B[i].length;
       B1[j].bits = Abits;

       if (Abits > B[i].bits)
       {
           freeBexp[j] = 1;
           B1[j].exps = (ulong *) flint_malloc(N*B[i].length*sizeof(ulong));
           mpoly_repack_monomials(B1[j].exps, Abits, B[i].exps, B[i].bits,
                                                           B[i].length, ctx->minfo);
       }
       else
       {
           freeBexp[j] = 0;
           B1[j].exps = B[i].exps;
       }
    }

    /* deal with aliasing and do multiplication */
    if (aliasing_required)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init3(T, maxlen, Abits, ctx);

        Alen = _fmpz_mpoly_addmul(T, B1, Blen, Abits, N, cmpmask, ctx);

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length_reset_bits(A, maxlen, Abits, ctx);

        Alen = _fmpz_mpoly_addmul(A, B1, Blen, Abits, N, cmpmask, ctx);
    }

    for (i = 0; i < Blen; i++)
    {
        if (freeBexp[i])
            flint_free(B1[i].exps);
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}



/* B is a list of Blen polynomials to be multiplied together and the result is placed in A.
 */

void fmpz_mpoly_addmul(
    fmpz_mpoly_t A,
    const fmpz_mpoly_struct * B,
    const slong Blen,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz * maxfields;
    fmpz * maxBfields;
    TMP_INIT;

    FLINT_ASSERT(Blen > 0);

    if (Blen == 1)
    {
        fmpz_mpoly_set(A, B + 0, ctx);
	return;
    }

    for (i = 0; i < Blen; i++)
    {
        if (B[i].length == 0)
        {
            fmpz_mpoly_zero(A, ctx);
            return;
        }
    }

    TMP_START;

    maxfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(maxfields + i);
        fmpz_init(maxBfields + i);
    }

    for (i = 0; i < Blen; i++)
    {
        mpoly_max_fields_fmpz(maxfields, (B+i)->exps, (B+i)->length, (B+i)->bits, ctx->minfo);
        _fmpz_vec_add(maxfields, maxfields, maxBfields, ctx->minfo->nfields);
    }

    _fmpz_mpoly_addmul_maxfields(A, B, Blen, maxBfields, ctx);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(maxfields + i);
        fmpz_clear(maxBfields + i);
    }

    TMP_END;
}

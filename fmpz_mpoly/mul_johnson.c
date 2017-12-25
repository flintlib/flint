/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product. This
   version of the function assumes the exponent vectors all fit in a
   single word. Assumes input polys are nonzero.
*/
slong _fmpz_mpoly_mul_johnson1(fmpz ** poly1, ulong ** exp1, slong * alloc,
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
slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                              mp_bitcnt_t bits, slong N, const ulong * cmpmask)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   slong * hind;
   int first, small;
   TMP_INIT;

   /* if exponent vectors fit in single word, call special version */
   if (N == 1)
      return _fmpz_mpoly_mul_johnson1(poly1, exp1, alloc,
                             poly2, exp2, len2, poly3, exp3, len3, cmpmask[0]);

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));
   for (i = 0; i < len2; i++)
      exp_list[i] = exps + i*N;

   /* space for heap indices */
   hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
   for (i = 0; i < len2; i++)
       hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   exp_next = 0;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

    if (bits <= FLINT_BITS)
        mpoly_monomial_add(heap[1].exp, exp2, exp3, N);
    else
        mpoly_monomial_add_mp(heap[1].exp, exp2, exp3, N);

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
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

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
               mpoly_monomial_set(e1 + k*N, exp, N);

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
               
               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

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

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
               exp_next--;
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

            if (bits <= FLINT_BITS)
                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);
            else
                mpoly_monomial_add_mp(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                      &next_loc, &heap_len, N, cmpmask))
               exp_next--;
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

void fmpz_mpoly_mul_johnson(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N, len = 0;
    mp_bitcnt_t exp_bits;
    fmpz * max_fields2, * max_fields3;
    ulong * cmpmask;
    ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
    int free2 = 0, free3 = 0;
    TMP_INIT;

    if (poly2->length == 0 || poly3->length == 0)
    {
        fmpz_mpoly_zero(poly1, ctx);
        return;
    }

    TMP_START;

    max_fields2 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    max_fields3 = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_init(max_fields2 + i);
        fmpz_init(max_fields3 + i);
    }
    mpoly_max_fields_fmpz(max_fields2, poly2->exps, poly2->length,
                                                      poly2->bits, ctx->minfo);
    mpoly_max_fields_fmpz(max_fields3, poly3->exps, poly3->length,
                                                      poly3->bits, ctx->minfo);
    _fmpz_vec_add(max_fields2, max_fields2, max_fields3, ctx->minfo->nfields);

    exp_bits = _fmpz_vec_max_bits(max_fields2, ctx->minfo->nfields);
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, exp_bits + 1);
    exp_bits = FLINT_MAX(exp_bits, poly2->bits);
    exp_bits = FLINT_MAX(exp_bits, poly3->bits);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nfields; i++)
    {
        fmpz_clear(max_fields2 + i);
        fmpz_clear(max_fields3 + i);
    }

    N = mpoly_words_per_exp(exp_bits, ctx->minfo);
    cmpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, exp_bits, ctx->minfo);

   /* ensure input exponents are packed into same sized fields as output */
   if (exp_bits > poly2->bits)
   {
      free2 = 1;
      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_repack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                    poly2->length, ctx->minfo);
   }

   if (exp_bits > poly3->bits)
   {
      free3 = 1;
      exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
      mpoly_repack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                    poly3->length, ctx->minfo);
   }

   /* deal with aliasing and do multiplication */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      /* algorithm more efficient if smaller poly first */
      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                                    exp_bits, N, cmpmask);
      else
         len = _fmpz_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                                    exp_bits, N, cmpmask);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      /* algorithm more efficient if smaller poly first */
      if (poly2->length > poly3->length)
         len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                                    exp_bits, N, cmpmask);
      else
         len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                                    exp_bits, N, cmpmask);
   }

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   _fmpz_mpoly_set_length(poly1, len, ctx);

   TMP_END;
}

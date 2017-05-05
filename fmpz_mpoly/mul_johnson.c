/*
    Copyright (C) 2016 William Hart

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

slong _fmpz_mpoly_mul_johnson1(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3)
{
   slong k;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong exp, cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, negate, small;
   TMP_INIT;

   TMP_START;

   small = _fmpz_mpoly_fits_small(poly2, len2) && _fmpz_mpoly_fits_small(poly3, len3);

   heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   
   next_free = 0;

   x = chain + next_free++;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0] + exp3[0], x);

   k = -WORD(1);

   while (heap_len > 1)
   {
      exp = heap[1].exp;
      k++;

      if (k >= *alloc)
      {
         p1 = (fmpz *) flint_realloc(p1, 2*sizeof(fmpz)*(*alloc));
         e1 = (ulong *) flint_realloc(e1, 2*sizeof(ulong)*(*alloc));
         flint_mpn_zero(p1 + *alloc, *alloc);
         (*alloc) *= 2;
      }

      first = 1;

      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && heap[1].exp == exp)
      {
         x = _mpoly_heap_pop1(heap, &heap_len);
         
         if (small)
         {
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));
               e1[k] = exp;
               first = 0; 
            } else
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

            }
      
            if (x->j < len3 - 1)
               Q[Q_len++] = x;

            while ((x = x->next) != NULL)
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               if (x->j < len3 - 1)
                  Q[Q_len++] = x;
            }
         } else
         {
            if (first)
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               e1[k] = exp;
               first = 0; 
            } else
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
      
            if (x->j < len3 - 1 || x->j == 0)
               Q[Q_len++] = x;

            while ((x = x->next) != NULL)
            {
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               if (x->j < len3 - 1 || x->j == 0)
                  Q[Q_len++] = x;
            }
         }
      }
      
      while (Q_len > 0)
      {
         x = Q[--Q_len];
     
         if (x->j == 0 && x->i < len2 - 1)
         {
            mpoly_heap_t * x2 = chain + next_free++;
            x2->i = x->i + 1;
            x2->j = 0;
            x2->next = NULL;

            _mpoly_heap_insert1(heap, exp2[x->i + 1] + exp3[0], x2, &heap_len);
         }

         if (x->j < len3 - 1)
         {
            x->j++;
            x->next = NULL;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &heap_len);
         }
      }     

      if (small)
      {
         negate = 0;

         if (0 > (slong) c[2])
         {
            c[0] = ~c[0];
            c[1] = ~c[1];
            c[2] = ~c[2];
            add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, 1);
            negate = 1;
         } 

         fmpz_set_ui(p1 + k, c[2]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[1]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[0]);
      
         if (negate)
            fmpz_neg(p1 + k, p1 + k);
      }
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                           const fmpz * poly3, const ulong * exp3, slong len3, slong N)
{
   slong i, k;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong * exp, * exps;
   ulong ** exp_list;
   slong exp_next;
   int first, negate, small;
   TMP_INIT;

   if (N == 1)
      return _fmpz_mpoly_mul_johnson1(poly1, exp1, alloc,
                                         poly2, exp2, len2, poly3, exp3, len3);

   TMP_START;

   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
   exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));

   for (i = 0; i < len2; i++)
      exp_list[i] = exps + i*N;

   next_free = 0;
   exp_next = 0;

   x = chain + next_free++;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

   mpoly_monomial_add(heap[1].exp, exp2, exp3, N);

   k = -WORD(1);

   while (heap_len > 1)
   {
      exp = heap[1].exp;
      k++;

      if (k >= *alloc)
      {
         p1 = (fmpz *) flint_realloc(p1, 2*sizeof(fmpz)*(*alloc));
         e1 = (ulong *) flint_realloc(e1, 2*N*sizeof(ulong)*(*alloc));
         flint_mpn_zero(p1 + *alloc, *alloc);
         (*alloc) *= 2;
      }

      first = 1;

      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N);
         
         if (small)
         {
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            if (x->j < len3 - 1)
               Q[Q_len++] = x;

            while ((x = x->next) != NULL)
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               if (x->j < len3 - 1)
                  Q[Q_len++] = x;
            }
         } else
         {
            if (first)
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
      
            if (x->j < len3 - 1)
               Q[Q_len++] = x;

            while ((x = x->next) != NULL)
            {
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               if (x->j < len3 - 1)
                  Q[Q_len++] = x;
            }
         }
      }
      
      while (Q_len > 0)
      {
         x = Q[--Q_len];
     
         if (x->j == 0 && x->i < len2 - 1)
         {
            mpoly_heap_t * x2 = chain + next_free++;
            x2->i = x->i + 1;
            x2->j = 0;
            x2->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp2 + (x->i + 1)*N, exp3, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len, N))
               exp_next--;
         }

         x->j++;
         x->next = NULL;

         mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

         if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
            exp_next--;
      }     
   
      if (small)
      {
         negate = 0;

         if (0 > (slong) c[2])
         {
            c[0] = ~c[0];
            c[1] = ~c[1];
            c[2] = ~c[2];
            add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, 1);
            negate = 1;
         } 

         fmpz_set_ui(p1 + k, c[2]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[1]);
         fmpz_mul_2exp(p1 + k, p1 + k, FLINT_BITS);
         fmpz_add_ui(p1 + k, p1 + k, c[0]);
      
         if (negate)
            fmpz_neg(p1 + k, p1 + k);
      }
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
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2;
   ulong * max_degs3;
   ulong max;
   ulong * exp2, * exp3;
   int free2, free3;

   TMP_INIT;

   if (poly2->length == 0 || poly3->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   max = 0;

   for (i = 0; i < ctx->n; i++)
   {
      max_degs3[i] += max_degs2[i];
      if (max_degs3[i] < max_degs2[i] || 0 > (slong) max_degs3[i]) 
         flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_johnson");

      if (max_degs3[i] > max)
         max = max_degs3[i];
   }

   bits = FLINT_BIT_COUNT(max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   exp3 = mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->length, ctx->n, poly3->bits);
   
   free3 = exp3 != poly3->exps;

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                           poly2->coeffs, exp2, poly2->length, N);
      else
         len = _fmpz_mpoly_mul_johnson(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length, N);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                           poly2->coeffs, exp2, poly2->length, N);
      else
         len = _fmpz_mpoly_mul_johnson(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length, N);
   }

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   _fmpz_mpoly_set_length(poly1, len, ctx);

   TMP_END;
}

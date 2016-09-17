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

int _fmpz_mpoly_fits_small(const fmpz * poly, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
   {
      if (COEFF_IS_MPZ(poly[i]))
         return 0;
   }

   return 1;
}

#define HEAP_LEFT(i) (2*(i))
#define HEAP_RIGHT(i) (2*(i) + 1)
#define HEAP_PARENT(i) ((i)/2)

#define HEAP_ASSIGN(h, c1, c2) \
   do {                                \
      (h).exp = (c1);                  \
      (h).next = (c2);                 \
   } while (0)

fmpz_mpoly_heap_t * _fmpz_mpoly_heap_pop(fmpz_mpoly_heap_s * heap, slong * heap_len)
{
   ulong exp;
   slong l, j, r, i = 1, len;
   fmpz_mpoly_heap_t * x = heap[1].next;

   (*heap_len)--;

   if ((len = *heap_len) > 1)
   {
      exp = heap[len].exp;

      while ((l = HEAP_LEFT(i)) <= len)
      {
         r = HEAP_RIGHT(i);
         j = r > len || heap[l].exp < heap[r].exp ? l : r;
         if (heap[j].exp < exp)
         {
            heap[i] = heap[j];
            i = j;
         } else
            break;
       }

       heap[i] = heap[len];
   }

   return x; 
}

void _fmpz_mpoly_heap_push(fmpz_mpoly_heap_s * heap, ulong exp, fmpz_mpoly_heap_t * x, slong * heap_len)
{
   slong j, i;
   
   i = *heap_len;
      
   (*heap_len)++;
   
   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (exp < heap[j].exp)
      {
         heap[i] = heap[j];
         i = j;
      } else
         break;
   }

   HEAP_ASSIGN(heap[i], exp, x);
}

void _fmpz_mpoly_heap_insert(fmpz_mpoly_heap_s * heap, ulong exp, fmpz_mpoly_heap_t * x, slong * heap_len)
{
   slong i = *heap_len, j, n = *heap_len;
   
   if (i != 1 && exp == heap[1].exp)
   {
      x->next = heap[1].next;
      heap[1].next = x;

      return;
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (exp == heap[j].exp)
      {
         x->next = heap[j].next;
         heap[j].next = x;

         return;
      } else if (exp < heap[j].exp)
         i = j;
      else
         break;
   }

   (*heap_len)++;

   while (n > i)
   {
      heap[n] = heap[HEAP_PARENT(n)];
      n = HEAP_PARENT(n);
   }

   HEAP_ASSIGN(heap[i], exp, x);
}

slong _fmpz_mpoly_mul_johnson1_si(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const slong * poly2, const ulong * exp2, slong len2,
                           const slong * poly3, const ulong * exp3, slong len3)
{
   slong i, k;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   fmpz_mpoly_heap_s * heap = (fmpz_mpoly_heap_s *) flint_malloc((len2 + 1)*sizeof(fmpz_mpoly_heap_s));
   fmpz_mpoly_heap_t * chain = (fmpz_mpoly_heap_t *) flint_malloc(len2*sizeof(fmpz_mpoly_heap_t));
   fmpz_mpoly_heap_t ** free_list = (fmpz_mpoly_heap_t **) flint_malloc(len2*sizeof(fmpz_mpoly_heap_t *));
   fmpz_mpoly_heap_t ** Q = (fmpz_mpoly_heap_t **) flint_malloc(len2*sizeof(fmpz_mpoly_heap_t *));
   fmpz_mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong exp, cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, negate;

   for (i = 0; i < len2; i++)
      free_list[i] = chain + i;

   next_free = 0;

   x = free_list[next_free++];
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
         x = _fmpz_mpoly_heap_pop(heap, &heap_len);
         
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
            c[2] = (p[1] >= 0) ? cy : (1 - cy);
         }
      
         if (x->j < len3 - 1)
            Q[Q_len++] = x;
         else
            free_list[--next_free] = x;

         while ((x = x->next) != NULL)
         {
            smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
            add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
            if (p[1] >= 0)
               c[2] += cy;
            else
               c[2] -= (1 - cy);

            if (x->j < len3 - 1)
               Q[Q_len++] = x;
            else
               free_list[--next_free] = x;
         }

         while (Q_len > 0)
         {
            x = Q[--Q_len];
     
            if (x->j == 0 && x->i < len2 - 1)
            {
               fmpz_mpoly_heap_t * x2 = free_list[next_free++];
               x2->i = x->i + 1;
               x2->j = 0;
               x2->next = NULL;

               _fmpz_mpoly_heap_push(heap, exp2[x->i + 1] + exp3[0], x2, &heap_len);
            }

            x->j++;
            x->next = NULL;
            _fmpz_mpoly_heap_insert(heap, exp2[x->i] + exp3[x->j], x, &heap_len);
         }     
      }

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

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   flint_free(heap);
   flint_free(chain);
   flint_free(free_list);
   flint_free(Q);

   return k;
}

slong _fmpz_mpoly_mul_johnson1(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                            const fmpz * poly3, const ulong * exp3, slong len3)
{
   if (_fmpz_mpoly_fits_small(poly2, len2) && _fmpz_mpoly_fits_small(poly3, len3))
      return _fmpz_mpoly_mul_johnson1_si(poly1, exp1, alloc,
                     (slong *) poly2, exp2, len2, (slong *) poly3, exp3, len3);
   else
      flint_throw(FLINT_ERROR, "Not implemented yet!");

   return 0;
}

void fmpz_mpoly_mul_johnson(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2;
   ulong * max_degs3;
   ulong max2 = 0, max3 = 0, max;
   ulong * exp2, * exp3;

   if (poly2->length == 0 || poly3->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return;
   }

   max_degs2 = fmpz_mpoly_max_degrees(poly2, ctx);
   max_degs3 = fmpz_mpoly_max_degrees(poly3, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max2)
         max2 = max_degs2[i];

      if (max_degs3[i] > max3)
         max3 = max_degs3[i];
   }

   max = max2 + max3;
   if (max < max2 || 0 > (slong) max)
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_johnson");

   bits = FLINT_BIT_COUNT(max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   exp2 = _fmpz_mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->bits, ctx->n, poly2->length);

   exp3 = _fmpz_mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->bits, ctx->n, poly3->length);

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   if (N == 1)
   {
      if (poly1 == poly2 || poly1 == poly3)
      {
         fmpz_mpoly_t temp;

         fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
         fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

         len = _fmpz_mpoly_mul_johnson1(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length);

         fmpz_mpoly_swap(temp, poly1, ctx);

         fmpz_mpoly_clear(temp, ctx);

      } else
      {
         fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
         fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

         len = _fmpz_mpoly_mul_johnson1(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                           poly3->coeffs, exp3, poly3->length);
      }
   } else
      flint_throw(FLINT_ERROR, "Not implemented yet!");

   if (exp2 != poly2->exps)
      flint_free(exp2);

   if (exp3 != poly3->exps)
      flint_free(exp3);

   flint_free(max_degs2);
   flint_free(max_degs3);

   _fmpz_mpoly_set_length(poly1, len, ctx);
}

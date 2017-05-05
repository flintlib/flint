/*
    Copyright (C) 2017 William Hart

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
#include "longlong.h"

slong _fmpz_mpoly_divides_monagan_pearce1(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
                const fmpz * poly3, const ulong * exp3, slong len3, slong bits)
{
   slong i, k, s;
   slong next_free, Q_len = 0, reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong exp;
   ulong c[3], p[2]; /* for accumulating coefficients */
   int first, d1, d2;
   ulong mask = 0, ub;
   fmpz_t mb, qc, r;
   int small;
   slong bits2, bits3;
   TMP_INIT;

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);
   fmpz_init(r);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
   chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
   Q = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   reuse = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));

   next_free = 0;

   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   k = -WORD(1);
   s = len3;
   
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[0], x);

   fmpz_neg(mb, poly3 + 0);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(poly3[0])/2 */
   
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
      d1 = 0;

      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && heap[1].exp == exp)
      {
         x = _mpoly_heap_pop1(heap, &heap_len);
         
         if (first)
         {
            d1 = mpoly_monomial_divides1(e1 + k, exp, exp3[0], mask);

            first = 0; 
         }

         if (small)
         {
            fmpz fc = poly2[x->j];

            if (x->i == -WORD(1))
            {
               if (!COEFF_IS_MPZ(fc))
               {
                  if (fc >= 0)
                     sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, fc);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, -fc);
               } else
               {
                  slong size = fmpz_size(poly2 + x->j);
                  __mpz_struct * m = COEFF_TO_PTR(fc);
                  if (fmpz_sgn(poly2 + x->j) < 0)
                     mpn_add(c, c, 3, m->_mp_d, size);
                  else
                     mpn_sub(c, c, 3, m->_mp_d, size);
               }
            } else
            {
               smul_ppmm(p[1], p[0], poly3[x->i], p1[x->j]);
               if (0 > (slong) p[1])
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
               else
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
            }
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + x->j);
            else
               fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
         }

         if (x->i != -WORD(1) || x->j < len2 - 1)
            Q[Q_len++] = x;
         else
            reuse[reuse_len++] = x;

         while ((x = x->next) != NULL)
         {
            if (small)
            {
               fmpz fc = poly2[x->j];

               if (x->i == -WORD(1))
               {
                  if (!COEFF_IS_MPZ(fc))
                  {
                     if (fc >= 0)
                        sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, fc);
                     else
                        add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, -fc);
                  } else
                  {
                     slong size = fmpz_size(poly2 + x->j);
                     __mpz_struct * m = COEFF_TO_PTR(fc);
                     if (fmpz_sgn(poly2 + x->j) < 0)
                        mpn_add(c, c, 3, m->_mp_d, size);
                     else
                        mpn_sub(c, c, 3, m->_mp_d, size);
                  }
               } else
               {
                  smul_ppmm(p[1], p[0], poly3[x->i], p1[x->j]);
                  if (0 > (slong) p[1])
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
               }
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + x->j);
               else
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
            }

            if (x->i != -WORD(1) || x->j < len2 - 1)
               Q[Q_len++] = x;
            else
               reuse[reuse_len++] = x;
         }
      }

      while (Q_len > 0)
      {
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            _mpoly_heap_insert1(heap, exp2[x->j], x, &heap_len);
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            _mpoly_heap_insert1(heap, exp3[x->i] + e1[x->j], x, &heap_len);
         } else if (x->j == k - 1)
         {
            s++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] == 0 && c[1] == 0 && c[0] == 0)) || (!small && fmpz_is_zero(qc)))
         k--;
      else
      {
         if (small)
         {
            ulong d[3];

            if (0 > (slong) c[2])
               mpn_neg(d, c, 3);
            else
               flint_mpn_copyi(d, c, 3);

            if (d[2] != 0 || ub < d[1] || (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
            {
               int negate = 0;

               /* get absolute value */
               if (0 > (slong) c[2])
               {
                  c[0] = ~c[0];
                  c[1] = ~c[1];
                  c[2] = ~c[2];
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, 1);
                  negate = 1;
               }

               /* set qc to abs(c) */
               fmpz_set_ui(qc, c[2]);
               fmpz_mul_2exp(qc, qc, FLINT_BITS);
               fmpz_add_ui(qc, qc, c[1]);
               fmpz_mul_2exp(qc, qc, FLINT_BITS);
               fmpz_add_ui(qc, qc, c[0]);

               /* correct sign */
               if (negate)
                  fmpz_neg(qc, qc);

               /* continue in non-small case from now on */
               fmpz_fdiv_qr(p1 + k, r, qc, mb);

               d2 = fmpz_is_zero(r);

               small = 0;
            } else /* quotient fits a small */
            {
               ulong r1;

               sdiv_qrnnd(p1[k], r1, c[1], c[0], *mb);
               
               d2 = r1 == 0;
            }
         } else /* quotient non-small case */
         {
            fmpz_fdiv_qr(p1 + k, r, qc, mb);

            d2 = fmpz_is_zero(r);
         }

         if (!d1 || !d2) /* not an exact division */
         {
            for (i = 0; i <= k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
         }

         for (i = 1; i < s; i++)
         {
            if (reuse_len != 0)
               x2 = reuse[--reuse_len];
            else
               x2 = chain + next_free++;
            
            x2->i = i;
            x2->j = k;
            x2->next = NULL;

            _mpoly_heap_insert1(heap, exp3[i] + e1[k], x2, &heap_len);
         }

         s = 1;
      } 
      
      fmpz_zero(qc);  
   }

   k++;

cleanup:

   fmpz_clear(mb);
   fmpz_clear(qc);
   fmpz_clear(r);

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1,
         slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2,
       const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N)
{
   slong i, k, s;
   slong next_free, Q_len = 0;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x, * x2;
   fmpz * p1 = *poly1;
   ulong * e1 = *exp1;
   ulong * exp, * exps;
   ulong ** exp_list;
   ulong c[3], p[2]; /* for accumulating coefficients */
   slong exp_next;
   int first, d1, d2;
   ulong mask = 0, ub;
   fmpz_t mb, qc, r;
   int small;
   slong bits2, bits3;
   TMP_INIT;

   if (N == 1)
      return _fmpz_mpoly_divides_monagan_pearce1(poly1, exp1, alloc,
                                    poly2, exp2, len2, poly3, exp3, len3, bits);

   TMP_START;

   fmpz_init(mb);
   fmpz_init(qc);
   fmpz_init(r);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   bits3 = _fmpz_vec_max_bits(poly3, len3);
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
   chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
   Q = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   reuse = (mpoly_heap_t **) TMP_ALLOC(len3*sizeof(mpoly_heap_t *));
   exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
   exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));

   for (i = 0; i < len3; i++)
      exp_list[i] = exps + i*N;

   next_free = 0;
   exp_next = 0;

   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   k = -WORD(1);
   s = len3;
   
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

   mpoly_monomial_set(heap[1].exp, exp2, N);

   fmpz_neg(mb, poly3 + 0);

   ub = ((ulong) FLINT_ABS(*mb)) >> 1; /* abs(poly3[0])/2 */
   
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
      d1 = 0;

      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N);
         
         if (first)
         {
            d1 = mpoly_monomial_divides(e1 + k*N, exp, exp3, N, mask);

            first = 0; 
         }

         if (small)
         {
            fmpz fc = poly2[x->j];

            if (x->i == -WORD(1))
            {
               if (!COEFF_IS_MPZ(fc))
               {
                  if (fc >= 0)
                     sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, fc);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, -fc);
               } else
               {
                  slong size = fmpz_size(poly2 + x->j);
                  __mpz_struct * m = COEFF_TO_PTR(fc);
                  if (fmpz_sgn(poly2 + x->j) < 0)
                     mpn_add(c, c, 3, m->_mp_d, size);
                  else
                     mpn_sub(c, c, 3, m->_mp_d, size);
               }
            } else
            {
               smul_ppmm(p[1], p[0], poly3[x->i], p1[x->j]);
               if (0 > (slong) p[1])
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
               else
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
            }
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + x->j);
            else
               fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
         }

         if (x->i != -WORD(1) || x->j < len2 - 1)
            Q[Q_len++] = x;
         else
            reuse[reuse_len++] = x;

         while ((x = x->next) != NULL)
         {
            if (small)
            {
               fmpz fc = poly2[x->j];

               if (x->i == -WORD(1))
               {
                  if (!COEFF_IS_MPZ(fc))
                  {
                     if (fc >= 0)
                        sub_dddmmmsss(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, fc);
                     else
                        add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, -fc);
                  } else
                  {
                     slong size = fmpz_size(poly2 + x->j);
                     __mpz_struct * m = COEFF_TO_PTR(fc);
                     if (fmpz_sgn(poly2 + x->j) < 0)
                        mpn_add(c, c, 3, m->_mp_d, size);
                     else
                        mpn_sub(c, c, 3, m->_mp_d, size);
                  }
               } else
               {
                  smul_ppmm(p[1], p[0], poly3[x->i], p1[x->j]);
                  if (0 > (slong) p[1])
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
               }
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + x->j);
               else
                  fmpz_addmul(qc, poly3 + x->i, p1 + x->j);
            }

            if (x->i != -WORD(1) || x->j < len2 - 1)
               Q[Q_len++] = x;
            else
               reuse[reuse_len++] = x;
         }
      }

      while (Q_len > 0)
      {
         x = Q[--Q_len];
     
         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_set(exp_list[exp_next], exp2 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
         } else if (x->j < k - 1)
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp3 + x->i*N, e1 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
         } else if (x->j == k - 1)
         {
            s++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] == 0 && c[1] == 0 && c[0] == 0)) || (!small && fmpz_is_zero(qc)))
         k--;
      else
      {
         if (small)
         {
            ulong d[3];

            if (0 > (slong) c[2])
               mpn_neg(d, c, 3);
            else
               flint_mpn_copyi(d, c, 3);

            if (d[2] != 0 || ub < d[1] || (ub == 0 && 0 > (slong) d[0])) /* quotient not a small */
            {
               int negate = 0;

               /* get absolute value */
               if (0 > (slong) c[2])
               {
                  c[0] = ~c[0];
                  c[1] = ~c[1];
                  c[2] = ~c[2];
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, 0, 1);
                  negate = 1;
               }

               /* set qc to abs(c) */
               fmpz_set_ui(qc, c[2]);
               fmpz_mul_2exp(qc, qc, FLINT_BITS);
               fmpz_add_ui(qc, qc, c[1]);
               fmpz_mul_2exp(qc, qc, FLINT_BITS);
               fmpz_add_ui(qc, qc, c[0]);

               /* correct sign */
               if (negate)
                  fmpz_neg(qc, qc);

               /* continue in non-small case from now on */
               fmpz_fdiv_qr(p1 + k, r, qc, mb);

               d2 = fmpz_is_zero(r);

               small = 0;
            } else /* quotient fits a small */
            {
               ulong r1;

               sdiv_qrnnd(p1[k], r1, c[1], c[0], *mb);
               
               d2 = r1 == 0;
            }
         } else /* quotient non-small case */
         {
            fmpz_fdiv_qr(p1 + k, r, qc, mb);

            d2 = fmpz_is_zero(r);
         }

         if (!d1 || !d2) /* not an exact division */
         {
            for (i = 0; i <= k; i++)
               _fmpz_demote(p1 + i);

            k = 0;

            goto cleanup;
         }

         for (i = 1; i < s; i++)
         {
            if (reuse_len != 0)
               x2 = reuse[--reuse_len];
            else
               x2 = chain + next_free++;
            
            x2->i = i;
            x2->j = k;
            x2->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp3 + i*N, e1 + k*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len, N))
               exp_next--;
         }

         s = 1;
      } 
      
      fmpz_zero(qc);  
   }

   k++;

cleanup:

   fmpz_clear(mb);
   fmpz_clear(qc);
   fmpz_clear(r);

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t poly1,
                  const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2, * max_degs3;
   ulong max = 0;
   ulong * exp2, * exp3, * expq;
   int free2 = 0, free3 = 0;
   ulong mask = 0;
   TMP_INIT;

   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divides_monagan_pearce");

   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return 1;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max)
         max = max_degs2[i];

      if (max_degs2[i] < max_degs3[i])
      {
         len = 0;

         goto cleanup;
      }
   }

   bits = FLINT_BIT_COUNT(max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   exp_bits = FLINT_MAX(exp_bits, poly2->bits);
   exp_bits = FLINT_MAX(exp_bits, poly3->bits);

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   expq = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   if (poly2->exps[poly2->length - 1] < poly3->exps[poly3->length - 1])
      goto cleanup;

   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   exp3 = mpoly_unpack_monomials(exp_bits, poly3->exps, 
                                           poly3->length, ctx->n, poly3->bits);
   
   free3 = exp3 != poly3->exps;

   for (i = 0; i < FLINT_BITS/exp_bits; i++)
      mask = (mask << exp_bits) + (UWORD(1) << (exp_bits - 1));

   /* check leading monomial divides */
   if (!mpoly_monomial_divides(expq, exp2 + (poly2->length - 1)*N,
                                        exp3 + (poly3->length - 1)*N, N, mask))
   {
      len = 0;

      goto cleanup;
   }

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      len = _fmpz_mpoly_divides_monagan_pearce(&temp->coeffs, &temp->exps,
                            &temp->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length/poly3->length + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      len = _fmpz_mpoly_divides_monagan_pearce(&poly1->coeffs, &poly1->exps,
                            &poly1->alloc, poly2->coeffs, exp2, poly2->length,
                              poly3->coeffs, exp3, poly3->length, exp_bits, N);
   }

cleanup:

   _fmpz_mpoly_set_length(poly1, len, ctx);

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   TMP_END;

   return (len != 0);
}

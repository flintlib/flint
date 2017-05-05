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

slong _fmpz_mpoly_divrem_ideal1(fmpz_mpoly_struct ** polyq, fmpz ** polyr,
  ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2,
 slong len2, const fmpz_mpoly_struct ** poly3, ulong * const * exp3, slong len,
                            slong bits, ulong maxn, const fmpz_mpoly_ctx_t ctx)
{
   slong i, l, n3;
   slong next_free, Q_len = 0, len3;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_nheap_t * chain;
   mpoly_nheap_t ** Q, ** reuse;
   mpoly_nheap_t * x, * x2;
   fmpz * p2 = *polyr;
   ulong * e2 = *expr;
   ulong exp, texp;
   ulong c[3], p[2]; /* for accumulating coefficients */
   ulong mask = 0;
   ulong * ub;
   slong * k, * s;
   fmpz_t qc, q, r;
   fmpz * mb;
   int small;
   slong bits2, bits3;
   int d1, d2, div_flag;
   TMP_INIT;

   TMP_START;

   fmpz_init(q);
   fmpz_init(qc);
   fmpz_init(r);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   
   bits3 = 0;
   len3 = 0;
   for (i = 0; i < len; i++)
   {
      slong b = fmpz_mpoly_max_bits(poly3[i]);
      bits3 = FLINT_MAX(bits3, b);
      len3 += poly3[i]->length;
   }
      
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap1_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap1_s));
   chain = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
   Q = (mpoly_nheap_t **) TMP_ALLOC(len3*sizeof(mpoly_nheap_t *));
   reuse = (mpoly_nheap_t **) TMP_ALLOC(len3*sizeof(mpoly_nheap_t *));
   k = (slong *) TMP_ALLOC(len*sizeof(slong));
   s = (slong *) TMP_ALLOC(len*sizeof(slong));
   ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
   mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

   next_free = 0;

   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   for (i = 0; i < len; i++)
   {
      k[i] = -WORD(1);
      s[i] = poly3[i]->length;
   }
   l = -WORD(1);
   
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->p = -WORD(1);
   x->next = NULL;

   HEAP_ASSIGN(heap[1], maxn - exp2[len2 - 1], x);

   for (i = 0; i < len; i++)
   {
      fmpz_init(mb + i);

      fmpz_neg(mb + i, poly3[i]->coeffs + poly3[i]->length - 1);

      ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
   }

   while (heap_len > 1)
   {
      exp = heap[1].exp;
      
      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && heap[1].exp == exp)
      {
         x = _mpoly_heap_pop1(heap, &heap_len);

         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;
               
         if (small)
         {
            fmpz fc = poly2[len2 - x->j - 1];

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
                  slong size = fmpz_size(poly2 + len2 - x->j - 1);
                  __mpz_struct * m = COEFF_TO_PTR(fc);
                  if (fmpz_sgn(poly2 + len2 - x->j - 1) < 0)
                     mpn_add(c, c, 3, m->_mp_d, size);
                  else
                     mpn_sub(c, c, 3, m->_mp_d, size);
               }
            } else
            {
               smul_ppmm(p[1], p[0], poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
               if (0 > (slong) p[1])
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
               else
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
            }
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
            else
               fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1, polyq[x->p]->coeffs + x->j);
         }

         if (x->i != -WORD(1) || x->j < len2 - 1)
            Q[Q_len++] = x;
         else
            reuse[reuse_len++] = x;

         while ((x = x->next) != NULL)
         {
            n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;

            if (small)
            {
               fmpz fc = poly2[len2 - x->j - 1];

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
                     slong size = fmpz_size(poly2 + len2 - x->j - 1);
                     __mpz_struct * m = COEFF_TO_PTR(fc);
                     if (fmpz_sgn(poly2 + len2 - x->j - 1) < 0)
                        mpn_add(c, c, 3, m->_mp_d, size);
                     else
                        mpn_sub(c, c, 3, m->_mp_d, size);
                  }
               } else
               {
                  smul_ppmm(p[1], p[0], poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
                  if (0 > (slong) p[1])
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
               }
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
               else
                  fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1, polyq[x->p]->coeffs + x->j);
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
         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;

         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            _mpoly_heap_insert1(heap, maxn - exp2[len2 - x->j - 1], x, &heap_len);
         } else if (x->j < k[x->p])
         {
            x->j++;
            x->next = NULL;

            _mpoly_heap_insert1(heap, maxn - exp3[x->p][n3 - x->i - 1] - polyq[x->p]->exps[x->j], x, &heap_len);
         } else if (x->j == k[x->p])
         {
            s[x->p]++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) || (!small && !fmpz_is_zero(qc)))
      {
         slong w;

         div_flag = 0;

         for (w = 0; w < len; w++)
         {
            n3 = poly3[w]->length;
            d1 = mpoly_monomial_divides1(&texp, maxn - exp3[w][n3 - 1], exp, mask);

            if (d1)
            {
               d2 = 0;

               if (small)
               {
                  ulong d[3];

                  if (0 > (slong) c[2])
                     mpn_neg(d, c, 3);
                  else
                     flint_mpn_copyi(d, c, 3);

                  if (d[2] != 0 || ub[w] < d[1] || (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
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

                     small = 0;
                  } else /* quotient fits a small */
                  {
                     ulong r1;
                     slong tq;

                     sdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);

                     d2 = (r1 == 0);

                     if (d2)
                     {
                        div_flag = 1;
                        k[w]++;

                        fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                        fmpz_set_si(polyq[w]->coeffs + k[w], tq);

                        polyq[w]->exps[k[w]] = texp;
                     }
                  }
               } 

               /* quotient non-small case */
               if (!small)
               {
                  fmpz_fdiv_qr(q, r, qc, mb + w);

                  d2 = fmpz_is_zero(r);
              
                  if (d2)
                  {
                     div_flag = 1;

                     k[w]++;

                     fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                     fmpz_set(polyq[w]->coeffs + k[w], q);
                     
                     polyq[w]->exps[k[w]] = texp;
                  }
               }

               if (d2)
               {
                  for (i = 1; i < s[w]; i++)
                  {
                     if (reuse_len != 0)
                        x2 = reuse[--reuse_len];
                     else
                        x2 = chain + next_free++;
            
                     x2->i = i;
                     x2->j = k[w];
                     x2->p = w;
                     x2->next = NULL;

                     _mpoly_heap_insert1(heap, maxn - exp3[w][n3 - i - 1] - polyq[w]->exps[k[w]], x2, &heap_len);
                  }
                  s[w] = 1;
                  break;
               }
            }
         }
         
         if (!div_flag)
         {
            l++;

            if (l >= *allocr)
            {
               p2 = (fmpz *) flint_realloc(p2, 2*sizeof(fmpz)*(*allocr));
               e2 = (ulong *) flint_realloc(e2, 2*sizeof(ulong)*(*allocr));
               flint_mpn_zero(p2 + *allocr, *allocr);
               (*allocr) *= 2;
            }

            if (small)
            {
               fmpz_set_si(p2 + l, c[0]);
               fmpz_neg(p2 + l, p2 + l);
            } else
               fmpz_neg(p2 + l, qc);

            e2[l] = maxn - exp;
         }
      } 
      
      fmpz_zero(qc);  
   }

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(r);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;
}

slong _fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** polyq, fmpz ** polyr,
  ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2,
 slong len2, const fmpz_mpoly_struct ** poly3, ulong * const * exp3, slong len,
                 slong N, slong bits, ulong * maxn, const fmpz_mpoly_ctx_t ctx)
{
   slong i, l, n3;
   slong next_free, Q_len = 0, len3;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_nheap_t * chain;
   mpoly_nheap_t ** Q, ** reuse;
   mpoly_nheap_t * x, * x2;
   fmpz * p2 = *polyr;
   ulong * e2 = *expr;
   ulong * exp, * exps, * texp;
   ulong ** exp_list;
   ulong c[3], p[2]; /* for accumulating coefficients */
   slong exp_next;
   ulong mask = 0;
   ulong * ub;
   slong * k, * s;
   fmpz_t qc, q, r;
   fmpz * mb;
   int small;
   slong bits2, bits3;
   int d1, d2, div_flag;
   TMP_INIT;

   TMP_START;

   fmpz_init(q);
   fmpz_init(qc);
   fmpz_init(r);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   
   bits3 = 0;
   len3 = 0;
   for (i = 0; i < len; i++)
   {
      slong b = fmpz_mpoly_max_bits(poly3[i]);
      bits3 = FLINT_MAX(bits3, b);
      len3 += poly3[i]->length;
   }
      
   /* allow one bit for sign, one bit for subtraction */
   small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) + FLINT_BITS - 2) &&
           FLINT_ABS(bits3) <= FLINT_BITS - 2;

   heap = (mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(mpoly_heap_s));
   chain = (mpoly_nheap_t *) TMP_ALLOC(len3*sizeof(mpoly_nheap_t));
   Q = (mpoly_nheap_t **) TMP_ALLOC(len3*sizeof(mpoly_nheap_t *));
   reuse = (mpoly_nheap_t **) TMP_ALLOC(len3*sizeof(mpoly_nheap_t *));
   exps = (ulong *) TMP_ALLOC(len3*N*sizeof(ulong));
   exp_list = (ulong **) TMP_ALLOC(len3*sizeof(ulong *));
   texp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   exp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   k = (slong *) TMP_ALLOC(len*sizeof(slong));
   s = (slong *) TMP_ALLOC(len*sizeof(slong));
   ub = (ulong *) TMP_ALLOC(len*sizeof(ulong));
   mb = (fmpz * ) TMP_ALLOC(len*sizeof(fmpz));

   for (i = 0; i < len3; i++)
      exp_list[i] = exps + i*N;

   next_free = 0;
   exp_next = 0;

   for (i = 0; i < FLINT_BITS/bits; i++)
      mask = (mask << bits) + (UWORD(1) << (bits - 1));

   for (i = 0; i < len; i++)
   {
      k[i] = -WORD(1);
      s[i] = poly3[i]->length;
   }
   l = -WORD(1);
   
   x = chain + next_free++;
   x->i = -WORD(1);
   x->j = 0;
   x->p = -WORD(1);
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

   mpoly_monomial_sub(heap[1].exp, maxn, exp2 + (len2 - 1)*N, N);

   for (i = 0; i < len; i++)
   {
      fmpz_init(mb + i);

      fmpz_neg(mb + i, poly3[i]->coeffs + poly3[i]->length - 1);

      ub[i] = ((ulong) FLINT_ABS(mb[i])) >> 1; /* abs(poly3[0])/2 */
   }

   while (heap_len > 1)
   {
      mpoly_monomial_set(exp, heap[1].exp, N);
      
      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N);

         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;
               
         if (small)
         {
            fmpz fc = poly2[len2 - x->j - 1];

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
                  slong size = fmpz_size(poly2 + len2 - x->j - 1);
                  __mpz_struct * m = COEFF_TO_PTR(fc);
                  if (fmpz_sgn(poly2 + len2 - x->j - 1) < 0)
                     mpn_add(c, c, 3, m->_mp_d, size);
                  else
                     mpn_sub(c, c, 3, m->_mp_d, size);
               }
            } else
            {
               smul_ppmm(p[1], p[0], poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
               if (0 > (slong) p[1])
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
               else
                  add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
            }
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
            else
               fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1, polyq[x->p]->coeffs + x->j);
         }

         if (x->i != -WORD(1) || x->j < len2 - 1)
            Q[Q_len++] = x;
         else
            reuse[reuse_len++] = x;

         while ((x = x->next) != NULL)
         {
            n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;

            if (small)
            {
               fmpz fc = poly2[len2 - x->j - 1];

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
                     slong size = fmpz_size(poly2 + len2 - x->j - 1);
                     __mpz_struct * m = COEFF_TO_PTR(fc);
                     if (fmpz_sgn(poly2 + len2 - x->j - 1) < 0)
                        mpn_add(c, c, 3, m->_mp_d, size);
                     else
                        mpn_sub(c, c, 3, m->_mp_d, size);
                  }
               } else
               {
                  smul_ppmm(p[1], p[0], poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
                  if (0 > (slong) p[1])
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], ~WORD(0), p[1], p[0]);
                  else
                     add_sssaaaaaa(c[2], c[1], c[0], c[2], c[1], c[0], 0, p[1], p[0]);
               }
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
               else
                  fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1, polyq[x->p]->coeffs + x->j);
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
         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;

         if (x->i == -WORD(1))
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_sub(exp_list[exp_next], maxn, exp2 + (len2 - x->j - 1)*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
         } else if (x->j < k[x->p])
         {
            x->j++;
            x->next = NULL;

            mpoly_monomial_sub(texp, maxn, exp3[x->p] + (n3 - x->i - 1)*N, N);
            mpoly_monomial_sub(exp_list[exp_next], texp, polyq[x->p]->exps + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
         } else if (x->j == k[x->p])
         {
            s[x->p]++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) || (!small && !fmpz_is_zero(qc)))
      {
         slong w;

         div_flag = 0;

         for (w = 0; w < len; w++)
         {
            n3 = poly3[w]->length;

            mpoly_monomial_sub(texp, maxn, exp3[w] + (n3 - 1)*N, N);

            d1 = mpoly_monomial_divides(texp, texp, exp, N, mask);

            if (d1)
            {
               d2 = 0;

               if (small)
               {
                  ulong d[3];

                  if (0 > (slong) c[2])
                     mpn_neg(d, c, 3);
                  else
                     flint_mpn_copyi(d, c, 3);

                  if (d[2] != 0 || ub[w] < d[1] || (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
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

                     small = 0;
                  } else /* quotient fits a small */
                  {
                     ulong r1;
                     slong tq;

                     sdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);

                     d2 = (r1 == 0);

                     if (d2)
                     {
                        div_flag = 1;
                        k[w]++;

                        fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                        fmpz_set_si(polyq[w]->coeffs + k[w], tq);

                        mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                     }
                  }
               } 

               /* quotient non-small case */
               if (!small)
               {
                  fmpz_fdiv_qr(q, r, qc, mb + w);

                  d2 = fmpz_is_zero(r);
              
                  if (d2)
                  {
                     div_flag = 1;

                     k[w]++;

                     fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                     fmpz_set(polyq[w]->coeffs + k[w], q);
                     
                     mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                  }
               }

               if (d2)
               {
                  for (i = 1; i < s[w]; i++)
                  {
                     if (reuse_len != 0)
                        x2 = reuse[--reuse_len];
                     else
                        x2 = chain + next_free++;
            
                     x2->i = i;
                     x2->j = k[w];
                     x2->p = w;
                     x2->next = NULL;

                     mpoly_monomial_sub(texp, maxn, exp3[w] + (n3 - i - 1)*N, N);
                     mpoly_monomial_sub(exp_list[exp_next], texp, polyq[w]->exps + k[w]*N, N);

                     if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len, N))
                        exp_next--;
                  }
                  s[w] = 1;
                  break;
               }
            }
         }
         
         if (!div_flag)
         {
            l++;

            if (l >= *allocr)
            {
               p2 = (fmpz *) flint_realloc(p2, 2*sizeof(fmpz)*(*allocr));
               e2 = (ulong *) flint_realloc(e2, 2*N*sizeof(ulong)*(*allocr));
               flint_mpn_zero(p2 + *allocr, *allocr);
               (*allocr) *= 2;
            }

            if (small)
            {
               fmpz_set_si(p2 + l, c[0]);
               fmpz_neg(p2 + l, p2 + l);
            } else
               fmpz_neg(p2 + l, qc);

            mpoly_monomial_sub(e2 + l*N, maxn, exp, N);
         }
      } 
      
      fmpz_zero(qc);  
   }

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(r);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;
}

void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
       const fmpz_mpoly_t poly2, const fmpz_mpoly_struct ** poly3, slong len,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, j, exp_bits, N, lenr = 0;
   slong len3 = 0;
   ulong * exp2;
   ulong ** exp3;
   ulong * max_degs3;
   int free2 = 0;
   int * free3;
   fmpz_mpoly_t temp2;
   fmpz_mpoly_struct * tr;
   ulong * maxn, * maxexp;
   int deg, rev;
   TMP_INIT;

   for (i = 0; i < len; i++)
   {  
      if (poly3[i]->length == 0)
         flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_ideal");

      len3 = FLINT_MAX(len3, poly3[i]->length);
   }

   if (poly2->length == 0)
   {
      for (i = 0; i < len; i++)
      {
         fmpz_mpoly_zero(q[i], ctx);
      }
      
      fmpz_mpoly_zero(r, ctx);

      return;
   }

   TMP_START;

   degrev_from_ord(deg, rev, ctx->ord);

   maxexp = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   free3 = (int *) TMP_ALLOC(len*sizeof(int));

   exp3 = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

   fmpz_mpoly_max_degrees(maxexp, poly2, ctx);
   
   exp_bits = poly2->bits;

   for (i = 0; i < len; i++)
   {
      fmpz_mpoly_max_degrees(max_degs3, poly3[i], ctx);

      for (j = 0; j < ctx->n; j++)
         maxexp[j] = FLINT_MAX(maxexp[j], max_degs3[j]);

      exp_bits = FLINT_MAX(exp_bits, poly3[i]->bits);
   }

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   maxn = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   mpoly_set_monomial(maxn, maxexp, exp_bits, ctx->n, deg, rev);

   exp2 = mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->length, ctx->n, poly2->bits);

   free2 = exp2 != poly2->exps;

   for (i = 0; i < len; i++)
   {
      exp3[i] = mpoly_unpack_monomials(exp_bits, poly3[i]->exps, 
                                     poly3[i]->length, ctx->n, poly3[i]->bits);
   
      free3[i] = exp3[i] != poly3[i]->exps;

      fmpz_mpoly_fit_length(q[i], 1, ctx);
      fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);
   }

   if (r == poly2)
   {
      fmpz_mpoly_init2(temp2, len3, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, len3, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);

      tr = r;
   }

   if (N == 1)
   {
      lenr = _fmpz_mpoly_divrem_ideal1(q, &tr->coeffs, &tr->exps,
                         &tr->alloc, poly2->coeffs, exp2, poly2->length,
                                       poly3, exp3, len, exp_bits, *maxn, ctx);
   } else
   {
      lenr = _fmpz_mpoly_divrem_ideal(q, &tr->coeffs, &tr->exps,
                         &tr->alloc, poly2->coeffs, exp2, poly2->length,
                                     poly3, exp3, len, N, exp_bits, maxn, ctx);
   }

   if (r == poly2)
   {
      fmpz_mpoly_swap(temp2, r, ctx);

      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(r, lenr, ctx);

   for (i = 0; i < len; i++)
      fmpz_mpoly_reverse(q[i], q[i], ctx);
   fmpz_mpoly_reverse(r, r, ctx);

   if (free2)
      flint_free(exp2);

   for (i = 0; i < len; i++)
   {
      if (free3[i])
         flint_free(exp3[i]);
   }

   TMP_END;
}

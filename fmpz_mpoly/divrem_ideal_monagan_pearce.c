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

/*
   As for divrem_monagan_pearce1 except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce1(fmpz_mpoly_struct ** polyq, fmpz ** polyr,
  ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2,
    slong len2, fmpz_mpoly_struct * const * poly3, ulong * const * exp3,
                 slong len, slong bits, ulong maxn, const fmpz_mpoly_ctx_t ctx)
{
   slong i, l, n3, w;
   slong next_free, Q_len = 0, len3;
   slong reuse_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_nheap_t * chain;
   mpoly_nheap_t ** Q, ** reuse;
   mpoly_nheap_t * x, * x2;
   fmpz * p2 = *polyr;
   ulong * e2 = *expr;
   ulong exp, texp;
   ulong c[3]; /* for accumulating coefficients */
   ulong mask = 0;
   ulong * ub;
   slong * k, * s;
   fmpz_t qc, q;
   fmpz * mb;
   int small;
   slong bits2, bits3;
   int d1, d2, div_flag;
   TMP_INIT;

   TMP_START;

   fmpz_init(q);
   fmpz_init(qc);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   
   bits3 = 0;
   len3 = 0;
   for (i = 0; i < len; i++)
   {
      slong b = FLINT_ABS(fmpz_mpoly_max_bits(poly3[i]));
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

      /* check there has been no overflow */
      if ((exp & mask) != 0)
      {
         for (i = 0; i < l; i++)
               _fmpz_demote(p2 + i);
         for (w = 0; w < len; w++)
         {
            for (i = 0; i < k[w]; i++)
               _fmpz_demote(polyq[w]->coeffs + i);

            k[w] = -WORD(1); /* we add 1 to this before exit */
         }

         l = -WORD(2); /* we add 1 to this upon return */

         goto cleanup;
      }

      c[0] = c[1] = c[2] = 0;

      while (heap_len > 1 && heap[1].exp == exp)
      {
         x = _mpoly_heap_pop1(heap, &heap_len);

         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;
               
         if (small)
         {
            if (x->i == -WORD(1))
               _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + len2 - x->j - 1);
            else
            {   _fmpz_mpoly_submul_uiuiui_fmpz(c,
                poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
            }
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
            else
               fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1,
                                                   polyq[x->p]->coeffs + x->j);
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
               if (x->i == -WORD(1))
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + len2 - x->j - 1);
               else
                  _fmpz_mpoly_submul_uiuiui_fmpz(c,
                poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
               else
                  fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1,
                                                   polyq[x->p]->coeffs + x->j);
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

            _mpoly_heap_insert1(heap, maxn - exp2[len2 - x->j - 1], x, &heap_len);
         } else if (x->j < k[x->p])
         {
            x->j++;
            x->next = NULL;

            n3 = poly3[x->p]->length;
            _mpoly_heap_insert1(heap, maxn - exp3[x->p][n3 - x->i - 1] -
                                        polyq[x->p]->exps[x->j], x, &heap_len);
         } else if (x->j == k[x->p])
         {
            s[x->p]++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
      {
         div_flag = 0;

         for (w = 0; w < len; w++)
         {
            n3 = poly3[w]->length;
            d1 = mpoly_monomial_divides1(&texp, maxn - 
                                                   exp3[w][n3 - 1], exp, mask);

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

                  if (d[2] != 0 || ub[w] <= d[1] ||
                   (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                  {
                     fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                     small = 0;
                  } else /* quotient fits a small */
                  {
                     slong r1;
                     slong tq;

                     sdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);

                     if (COEFF_IS_MPZ(FLINT_ABS(tq))) /* quotient too large */
                     {
                        /* upgrade to multiprecision accumulated coeffs */
                        fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                        small = 0;
                     } else
                     {
                        div_flag = (r1 == 0);

                        d2 = tq != 0;

                        if (d2)
                        {
                           k[w]++;

                           fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                           fmpz_set_si(polyq[w]->coeffs + k[w], tq);

                           polyq[w]->exps[k[w]] = texp;
                        }

                        c[0] = r1;
                        c[2] = c[1] = r1 < 0 ? ~WORD(0) : WORD(0);
                     }
                  }
               } 

               /* quotient non-small case */
               if (!small)
               {
                  fmpz_fdiv_qr(q, qc, qc, mb + w);

                  div_flag = fmpz_is_zero(qc);

                  d2 = !fmpz_is_zero(q);

                  if (d2)
                  {
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

                     _mpoly_heap_insert1(heap, maxn - exp3[w][n3 - i - 1] -
                                          polyq[w]->exps[k[w]], x2, &heap_len);
                  }
                  s[w] = 1;
               }
            }
         }
         
         if (!div_flag)
         {
            l++;

            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, 1);

            if (small)
            {
               fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
               fmpz_neg(p2 + l, p2 + l);
            } else
               fmpz_neg(p2 + l, qc);

            e2[l] = maxn - exp;
         }
      } 
      
      fmpz_zero(qc);  
   }

cleanup:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;
}

/*
   As for divrem_monagan_pearce except that an array of divisor polynomials is
   passed and an array of quotient polynomials is returned. These are not in
   low level format.
*/
slong _fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq, fmpz ** polyr,
  ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2,
    slong len2, fmpz_mpoly_struct * const * poly3, ulong * const * exp3, 
      slong len, slong N, slong bits, ulong * maxn, const fmpz_mpoly_ctx_t ctx)
{
   slong i, l, n3, w;
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
   ulong c[3]; /* for accumulating coefficients */
   slong exp_next;
   ulong mask = 0;
   ulong * ub;
   slong * k, * s;
   fmpz_t qc, q;
   fmpz * mb;
   int small;
   slong bits2, bits3;
   int d1, d2, div_flag;
   TMP_INIT;

   if (N == 1)
      return _fmpz_mpoly_divrem_ideal_monagan_pearce1(polyq, polyr, expr, allocr,
                        poly2, exp2, len2, poly3, exp3, len, bits, *maxn, ctx);

   TMP_START;

   fmpz_init(q);
   fmpz_init(qc);

   bits2 = _fmpz_vec_max_bits(poly2, len2);
   
   bits3 = 0;
   len3 = 0;
   for (i = 0; i < len; i++)
   {
      slong b = FLINT_ABS(fmpz_mpoly_max_bits(poly3[i]));
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

      /* check there has been no overflow */
      if (mpoly_monomial_overflows(exp, N, mask))
      {
         for (i = 0; i < l; i++)
               _fmpz_demote(p2 + i);
         for (w = 0; w < len; w++)
         {
            for (i = 0; i < k[w]; i++)
               _fmpz_demote(polyq[w]->coeffs + i);

            k[w] = -WORD(1); /* we add 1 to this before exit */
         }

         l = -WORD(2); /* we add 1 to this upon return */

         goto cleanup2;
      }

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N);

         n3 = x->i == -WORD(1) ? 0 : poly3[x->p]->length;
               
         if (small)
         {
            if (x->i == -WORD(1))
               _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + len2 - x->j - 1);
            else
               _fmpz_mpoly_submul_uiuiui_fmpz(c,
                poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
         } else
         {
            if (x->i == -WORD(1))
               fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
            else
               fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1,
                                                   polyq[x->p]->coeffs + x->j);
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
               if (x->i == -WORD(1))
                  _fmpz_mpoly_sub_uiuiui_fmpz(c, poly2 + len2 - x->j - 1);
               else
                  _fmpz_mpoly_submul_uiuiui_fmpz(c,
                poly3[x->p]->coeffs[n3 - x->i - 1], polyq[x->p]->coeffs[x->j]);
            } else
            {
               if (x->i == -WORD(1))
                  fmpz_sub(qc, qc, poly2 + len2 - x->j - 1);
               else
                  fmpz_addmul(qc, poly3[x->p]->coeffs + n3 - x->i - 1,
                                                   polyq[x->p]->coeffs + x->j);
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

            mpoly_monomial_sub(exp_list[exp_next], maxn, 
                                                exp2 + (len2 - x->j - 1)*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                                                 &heap_len, N))
               exp_next--;
         } else if (x->j < k[x->p])
         {
            x->j++;
            x->next = NULL;

            n3 = poly3[x->p]->length;
            mpoly_monomial_sub(texp, maxn, exp3[x->p] + (n3 - x->i - 1)*N, N);
            mpoly_monomial_sub(exp_list[exp_next],
                                          texp, polyq[x->p]->exps + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, 
                                                                 &heap_len, N))
               exp_next--;
         } else if (x->j == k[x->p])
         {
            s[x->p]++;
            reuse[reuse_len++] = x;
         }
      }

      if ((small && (c[2] != 0 || c[1] != 0 || c[0] != 0)) ||
                                                 (!small && !fmpz_is_zero(qc)))
      {
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

                  if (d[2] != 0 || ub[w] <= d[1] ||
                   (ub[w] == 0 && 0 > (slong) d[0])) /* quotient not a small */
                  {
                     fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                     small = 0;
                  } else /* quotient fits a small */
                  {
                     slong r1;
                     slong tq;

                     sdiv_qrnnd(tq, r1, c[1], c[0], mb[w]);

                     if (COEFF_IS_MPZ(FLINT_ABS(tq))) /* quotient too large */
                     {
                        /* upgrade to multiprecision accumulated coeffs */
                        fmpz_set_signed_uiuiui(qc, c[2], c[1], c[0]);

                        small = 0;
                     } else
                     {
                        d2 = tq != 0;

                        div_flag = r1 == 0;

                        if (d2)
                        {
                           k[w]++;

                           fmpz_mpoly_fit_length(polyq[w], k[w] + 1, ctx);

                           fmpz_set_si(polyq[w]->coeffs + k[w], tq);

                           mpoly_monomial_set(polyq[w]->exps + k[w]*N, texp, N);
                        }

                        c[0] = r1;
                        c[2] = c[1] = r1 < 0 ? ~WORD(0) : WORD(0);
                     }
                  }
               } 

               /* quotient non-small case */
               if (!small)
               {
                  fmpz_fdiv_qr(q, qc, qc, mb + w);

                  d2 = !fmpz_is_zero(q);

                  div_flag = fmpz_is_zero(qc);

                  if (d2)
                  {
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

                     mpoly_monomial_sub(texp, maxn, exp3[w] + 
                                                            (n3 - i - 1)*N, N);
                     mpoly_monomial_sub(exp_list[exp_next], texp, 
                                                   polyq[w]->exps + k[w]*N, N);

                     if (!_mpoly_heap_insert(heap, exp_list[exp_next++],
                                                             x2, &heap_len, N))
                        exp_next--;
                  }
                  s[w] = 1;
               }
            }
         }
         
         if (!div_flag)
         {
            l++;

            _fmpz_mpoly_fit_length(&p2, &e2, allocr, l + 1, N);

            if (small)
            {
               fmpz_set_signed_uiuiui(p2 + l, c[2], c[1], c[0]);
               fmpz_neg(p2 + l, p2 + l);
            } else
               fmpz_neg(p2 + l, qc);

            mpoly_monomial_sub(e2 + l*N, maxn, exp, N);
         }
      } 
      
      fmpz_zero(qc);  
   }

cleanup2:

   for (i = 0; i < len; i++)
      _fmpz_mpoly_set_length(polyq[i], k[i] + 1, ctx); 

   for (i = 0; i < len; i++)
      fmpz_clear(mb + i);
   fmpz_clear(qc);
   fmpz_clear(q);

   (*polyr) = p2;
   (*expr) = e2;
   
   TMP_END;

   return l + 1;
}

/* Assumes divisor polys don't alias any output polys */
void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r,
    const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong i, j, exp_bits, N, lenr = 0;
   slong len3 = 0;
   ulong * exp2;
   ulong ** exp3;
   int free2 = 0;
   int * free3;
   fmpz_mpoly_t temp2;
   fmpz_mpoly_struct * tr;
   ulong * maxn, * maxexp;
   TMP_INIT;

   /* check none of the divisor polynomials is zero */
   for (i = 0; i < len; i++)
   {  
      if (poly3[i]->length == 0)
         flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_divrem_ideal_monagan_pearce");

      len3 = FLINT_MAX(len3, poly3[i]->length);
   }

   /* dividend is zero, write out quotients and remainder */
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

   maxexp = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   free3 = (int *) TMP_ALLOC(len*sizeof(int));

   exp3 = (ulong **) TMP_ALLOC(len*sizeof(ulong *));

   exp_bits = poly2->bits;

   /* compute maximum degrees that can occur in any input or output polys */
   for (i = 0; i < len; i++)
      exp_bits = FLINT_MAX(exp_bits, poly3[i]->bits);

   /* compute bound for degrees appearing in inputs and outputs for each var */
   for (i = 0; i < ctx->n; i++)
      maxexp[i] = (UWORD(1) << (exp_bits - 1)) - UWORD(1);

   /* number of words in exponent vector */
   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   /* pack maxexp monomial vector */
   maxn = (ulong *) TMP_ALLOC(8*N*sizeof(ulong));

   mpoly_set_monomial(maxn, maxexp, exp_bits, ctx->n, 0, 0);

   /* ensure input exponents packed to same size as output exponents */
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
      q[i]->bits = exp_bits;
   }

   /* check leading mon. of at least one divisor is at most that of dividend */
   for (i = 0; i < len; i++)
   {
      if (!mpoly_monomial_lt(exp3[i] + (poly3[i]->length - 1)*N, 
                         exp2 + (poly2->length - 1)*N, N))
         break;
   }

   if (i == len)
   {
      fmpz_mpoly_set(r, poly2, ctx);
      for (i = 0; i < len; i++)
         fmpz_mpoly_zero(q[i], ctx);

      goto cleanup3;
   }

   /* take care of aliasing */
   if (r == poly2)
   {
      fmpz_mpoly_init2(temp2, len3, ctx);
      fmpz_mpoly_fit_bits(temp2, exp_bits, ctx);
      temp2->bits = exp_bits;

      tr = temp2;
   } else
   {
      fmpz_mpoly_fit_length(r, len3, ctx);
      fmpz_mpoly_fit_bits(r, exp_bits, ctx);
      r->bits = exp_bits;

      tr = r;
   }

   /* do division with remainder */
   while (exp_bits <= FLINT_BITS)
   {
      ulong * old_exp2 = exp2, * old_exp3;

      lenr = _fmpz_mpoly_divrem_ideal_monagan_pearce(q, &tr->coeffs, &tr->exps,
                         &tr->alloc, poly2->coeffs, exp2, poly2->length,
                                     poly3, exp3, len, N, exp_bits, maxn, ctx);

      if (lenr >= 0) /* check if division was successful */
         break;

      exp_bits *= 2;

      if (exp_bits > FLINT_BITS)
         break;

      exp2 = mpoly_unpack_monomials(exp_bits, old_exp2, 
                                            poly2->length, ctx->n, exp_bits/2);

      if (free2)
         flint_free(old_exp2);

      free2 = 1;
 
      fmpz_mpoly_fit_bits(tr, exp_bits, ctx);

      tr->bits = exp_bits;

      for (i = 0; i < len; i++)
      {
         old_exp3 = exp3[i];

         exp3[i] = mpoly_unpack_monomials(exp_bits, old_exp3, 
                                         poly3[i]->length, ctx->n, exp_bits/2);
   
         if (free3[i])
            flint_free(old_exp3);

         free3[i] = 1; 

         for (j = 0; j < ctx->n; j++)
            maxexp[i] = (UWORD(1) << (exp_bits - 1)) - UWORD(1);

         fmpz_mpoly_fit_bits(q[i], exp_bits, ctx);

         q[i]->bits = exp_bits;
      }

      mpoly_set_monomial(maxn, maxexp, exp_bits, ctx->n, 0, 0);

      N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;
   }

   if (lenr < 0)
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_divrem_ideal_monagan_pearce");

   /* take care of aliasing */
   if (r == poly2)
   {
      fmpz_mpoly_swap(temp2, r, ctx);

      fmpz_mpoly_clear(temp2, ctx);
   } 

   _fmpz_mpoly_set_length(r, lenr, ctx);

   /* quotient polys and remainder poly are reversed, so correct order */
   for (i = 0; i < len; i++)
      fmpz_mpoly_reverse(q[i], q[i], ctx);
   fmpz_mpoly_reverse(r, r, ctx);

cleanup3:

   if (free2)
      flint_free(exp2);

   for (i = 0; i < len; i++)
   {
      if (free3[i])
         flint_free(exp3[i]);
   }

   TMP_END;
}

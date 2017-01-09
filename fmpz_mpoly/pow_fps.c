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

slong _fmpz_mpoly_pow_fps1(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2, slong k)
{
   const slong topbit = (WORD(1) << (FLINT_BITS - 1));
   const slong mask = ~topbit;
   slong i, rnext, g_alloc, gnext;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap1_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1, * gc = NULL;
   ulong * e1 = *exp1, * ge, * fik;
   ulong exp, finalexp, temp2;
   slong * largest;
   fmpz_t t1, C, S, temp1;
   int first;
   TMP_INIT;

   TMP_START;

   heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
   /* 2x as we pull from heap and insert more before processing pulled ones */
   chain = (mpoly_heap_t *) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t));
   reuse = (mpoly_heap_t **) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t *));
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   /* we add 1 to all entries of largest to free up the value 0 */
   largest = (slong *) TMP_ALLOC(len2*sizeof(slong));

   fmpz_init(t1);
   fmpz_init(C);
   fmpz_init(S);
   fmpz_init(temp1);

   for (i = 0; i < 2*len2; i++)
      reuse[i] = chain + i;

   g_alloc = k*(len2 - 1) + 1;
   ge = (ulong *) flint_malloc(g_alloc*sizeof(ulong));
   gnext = 0;
   rnext = 0;

   ge[0] = exp2[0]*(k - 1);

   e1[0] = exp2[0]*k;

   gc = (fmpz *) flint_calloc(g_alloc, sizeof(fmpz));
   fmpz_pow_ui(gc + 0, poly2 + 0, k - 1);
   fmpz_mul(p1 + 0, gc + 0, poly2 + 0);

   next_free = 0;

   x = reuse[next_free++];
   x->i = 1;
   x->j = 0;
   x->next = NULL;

   HEAP_ASSIGN(heap[1], exp2[1] + ge[0], x);

   for (i = 0; i < len2; i++)
      largest[i] = topbit;
   largest[1] = 1;

   fik = (ulong *) TMP_ALLOC(len2*sizeof(ulong));

   for (i = 0; i < len2; i++)
      fik[i] = exp2[i]*(k - 1);

   finalexp = exp2[len2 - 1]*(k - 1) + exp2[0];

   while (heap_len > 1)
   {
      exp = heap[1].exp;

      rnext++;
      gnext++;

      if (rnext >= *alloc)
      {
         p1 = (fmpz *) flint_realloc(p1, 2*sizeof(fmpz)*(*alloc));
         e1 = (ulong *) flint_realloc(e1, 2*sizeof(ulong)*(*alloc));
         flint_mpn_zero(p1 + *alloc, *alloc);
         (*alloc) *= 2;
      }

      if (gnext >= g_alloc)
      {
         ge = (ulong *) flint_realloc(ge, 2*sizeof(ulong)*g_alloc);
         gc = (fmpz *) flint_realloc(gc, 2*sizeof(fmpz)*g_alloc);
         flint_mpn_zero(gc + g_alloc, g_alloc);
         g_alloc *= 2;
      }

      first = 1;

      fmpz_zero(C);
      fmpz_zero(S);

      while (heap_len > 1 && heap[1].exp == exp)
      {
         x = _mpoly_heap_pop1(heap, &heap_len);

         largest[x->i] |= topbit;

         fmpz_mul(t1, poly2 + x->i, gc + x->j);
         fmpz_add(S, S, t1);

         if (exp <= finalexp)
         {
            temp2 = fik[x->i] - ge[x->j];
            if ((slong) temp2 < 0)
               fmpz_submul_ui(C, t1, -temp2);
            else
               fmpz_addmul_ui(C, t1, temp2);
         }

         if (first)
         {
            ge[gnext] = exp - exp2[0];
            first = 0; 
         }
      
         Q[Q_len++] = x;

         while ((x = x->next) != NULL)
         {
            largest[x->i] |= topbit;

            fmpz_mul(t1, poly2 + x->i, gc + x->j);
            fmpz_add(S, S, t1);

            if (exp <= finalexp)
            {
               temp2 = fik[x->i] - ge[x->j];
               if ((slong) temp2 < 0)
                  fmpz_submul_ui(C, t1, -temp2);
               else
                  fmpz_addmul_ui(C, t1, temp2);
            }

            Q[Q_len++] = x;
         }
      }
      
      while (Q_len > 0)
      {
         slong i, j;

         x = Q[--Q_len];
         i = x->i;
         j = x->j;

         if (i < len2 - 1 && largest[i + 1] == (j | topbit))
         {
            x->i++;
            x->next = NULL;

            _mpoly_heap_insert1(heap, exp2[i + 1] + ge[j], x, &heap_len);
            largest[i + 1] = j + 1;
         } else
            reuse[--next_free] = x;

         if (j < gnext - 1 && (largest[i] & mask) < j + 2)
         {
            x = reuse[next_free++];

            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            _mpoly_heap_insert1(heap, exp2[i] + ge[j + 1], x, &heap_len);
            largest[i] = j + 2;
         }
      }

      if (!fmpz_is_zero(C))
      {
         fmpz_divexact_ui(temp1, C, exp - k*exp2[0]);
         fmpz_add(S, S, temp1);
         fmpz_divexact(gc + gnext, temp1, poly2 + 0);

         if ((largest[1] & topbit) != 0)
         {
            x = reuse[next_free++];

            x->i = 1;
            x->j = gnext;
            x->next = NULL;

            _mpoly_heap_insert1(heap, exp2[1] + ge[gnext], x, &heap_len);
            largest[1] = gnext + 1;
         }
      }

      if (!fmpz_is_zero(S))
      {
         fmpz_set(p1 + rnext, S);
         e1[rnext] = ge[gnext] + exp2[0];
      } else
         rnext--;

      if (fmpz_is_zero(C))
         gnext--;
   }

   rnext++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   fmpz_clear(t1);
   fmpz_clear(C);
   fmpz_clear(S);
   fmpz_clear(temp1);

   flint_free(ge);
   for (i = 0; i < g_alloc; i++)
      fmpz_clear(gc + i);
   flint_free(gc);

   TMP_END;

   return rnext;
}

void fmpz_set_mpn(fmpz_t f, ulong * c_in, slong n)		
{		
    slong i;
    ulong * c;
    
    TMP_INIT;

    TMP_START;

    c = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    for (i = 0; i < n; i++)
       c[i] = c_in[n - i - 1];

    while (n > 0 && c[n - 1] == 0)		
       n--;		
 		
    if (n <= 1)		
       fmpz_set_ui(f, c[0]);		
    else		
    {		
       __mpz_struct * mpz = _fmpz_promote(f);		
 	
       mpz_realloc2(mpz, n*FLINT_BITS);		
 		
       mpn_copyi(mpz->_mp_d, c, n);		
       mpz->_mp_size = n;		
    }	

    TMP_END;
 }		
 

void fmpz_set_mpn_signed(fmpz_t f, ulong * c_in, slong n)		
{		
    slong i;
    ulong * c;
    int neg;

    TMP_INIT;

    TMP_START;

    c = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    for (i = 0; i < n; i++)
       c[i] = c_in[n - i - 1];

    neg = 0 > (slong) c[n - 1];		
		
    if (neg)		
       mpn_neg_n(c, c, n);		
		
    while (n > 0 && c[n - 1] == 0)		
       n--;		
 		
    if (n <= 1)		
       fmpz_set_ui(f, c[0]);		
    else		
    {		
       __mpz_struct * mpz = _fmpz_promote(f);		
 	
       mpz_realloc2(mpz, n*FLINT_BITS);		
 		
       mpn_copyi(mpz->_mp_d, c, n);		
       mpz->_mp_size = n;		
    }		
 		
    if (neg)		
       fmpz_neg(f, f);

    TMP_END;		
 }		
 

slong _fmpz_mpoly_pow_fps(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2, slong k, slong N)
{
   const slong topbit = (WORD(1) << (FLINT_BITS - 1));
   const slong mask = ~topbit;
   slong i, rnext, g_alloc, gnext, exp_next;
   slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
   mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   mpoly_heap_t ** Q, ** reuse;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1, * gc = NULL;
   ulong * e1 = *exp1, * ge, * fik, * exp, * exps;
   ulong ** exp_list;
   ulong * finalexp, * temp2;
   slong * largest;
   fmpz_t t1, t2, C, S, temp1;
   int first;
   TMP_INIT;

   TMP_START;

   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* 2x as we pull from heap and insert more before processing pulled ones */
   chain = (mpoly_heap_t *) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t));
   reuse = (mpoly_heap_t **) TMP_ALLOC(2*len2*sizeof(mpoly_heap_t *));
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   /* we add 1 to all entries of largest to free up the value 0 */
   largest = (slong *) TMP_ALLOC(len2*sizeof(slong));
   exps = (ulong *) TMP_ALLOC((len2 + 1)*N*sizeof(ulong));
   exp_list = (ulong **) TMP_ALLOC((len2 + 1)*sizeof(ulong *));
   finalexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
   temp2 = (ulong *) TMP_ALLOC(N*sizeof(ulong));

   fmpz_init(t1);
   fmpz_init(t2);
   fmpz_init(C);
   fmpz_init(S);
   fmpz_init(temp1);

   for (i = 0; i < len2 + 1; i++)
      exp_list[i] = exps + i*N;

   exp_next = 0;

   for (i = 0; i < 2*len2; i++)
      reuse[i] = chain + i;

   g_alloc = k*(len2 - 1) + 1;
   ge = (ulong *) flint_malloc(g_alloc*sizeof(ulong)*N);
   gnext = 0;
   rnext = 0;

   mpoly_monomial_mul_si(ge + 0, exp2 + 0, N, k - 1);

   mpoly_monomial_mul_si(e1 + 0, exp2 + 0, N, k);

   gc = (fmpz *) flint_calloc(g_alloc, sizeof(fmpz));
   fmpz_pow_ui(gc + 0, poly2 + 0, k - 1);
   fmpz_mul(p1 + 0, gc + 0, poly2 + 0);

   next_free = 0;

   x = reuse[next_free++];
   x->i = 1;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];
   
   mpoly_monomial_add(heap[1].exp, exp2 + N, ge + 0, N);

   for (i = 0; i < len2; i++)
      largest[i] = topbit;
   largest[1] = 1;

   fik = (ulong *) TMP_ALLOC(N*len2*sizeof(ulong));

   for (i = 0; i < len2; i++)
      mpoly_monomial_mul_si(fik + i*N, exp2 + i*N, N, k - 1);

   mpoly_monomial_mul_si(finalexp, exp2 + (len2 - 1)*N, N, k - 1);
   mpoly_monomial_add(finalexp, finalexp, exp2 + 0, N);

   while (heap_len > 1)
   {
      exp = heap[1].exp;

      rnext++;
      gnext++;

      if (rnext >= *alloc)
      {
         p1 = (fmpz *) flint_realloc(p1, 2*sizeof(fmpz)*(*alloc));
         e1 = (ulong *) flint_realloc(e1, 2*N*sizeof(ulong)*(*alloc));
         flint_mpn_zero(p1 + *alloc, *alloc);
         (*alloc) *= 2;
      }

      if (gnext >= g_alloc)
      {
         ge = (ulong *) flint_realloc(ge, 2*N*sizeof(ulong)*g_alloc);
         gc = (fmpz *) flint_realloc(gc, 2*sizeof(fmpz)*g_alloc);
         flint_mpn_zero(gc + g_alloc, g_alloc);
         g_alloc *= 2;
      }

      first = 1;

      fmpz_zero(C);
      fmpz_zero(S);

      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         x = _mpoly_heap_pop(heap, &heap_len, N);

         largest[x->i] |= topbit;

         fmpz_mul(t1, poly2 + x->i, gc + x->j);
         fmpz_add(S, S, t1);

         if (!mpoly_monomial_lt(finalexp, exp, N))
         {
            mpoly_monomial_sub(temp2, fik + x->i*N, ge + x->j*N, N);
            fmpz_set_mpn_signed(t2, temp2, N);

            fmpz_addmul(C, t1, t2);
         }

         if (first)
         {
            mpoly_monomial_sub_no_borrow(ge + gnext*N, exp, exp2 + 0, N);
            first = 0; 
         }
      
         Q[Q_len++] = x;

         while ((x = x->next) != NULL)
         {
            largest[x->i] |= topbit;

            fmpz_mul(t1, poly2 + x->i, gc + x->j);
            fmpz_add(S, S, t1);

            if (!mpoly_monomial_lt(finalexp, exp, N))
            {
               mpoly_monomial_sub(temp2, fik + x->i*N, ge + x->j*N, N);
               fmpz_set_mpn_signed(t2, temp2, N);

               fmpz_addmul(C, t1, t2);
            }

            Q[Q_len++] = x;
         }
      }
      
      while (Q_len > 0)
      {
         slong i, j;

         x = Q[--Q_len];
         i = x->i;
         j = x->j;

         if (i < len2 - 1 && largest[i + 1] == (j | topbit))
         {
            x->i++;
            x->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp2 + (i + 1)*N, ge + j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
 
            largest[i + 1] = j + 1;
         } else
            reuse[--next_free] = x;

         if (j < gnext - 1 && (largest[i] & mask) < j + 2)
         {
            x = reuse[next_free++];

            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp2 + i*N, ge + (j + 1)*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
 
            largest[i] = j + 2;
         }
      }

      if (!fmpz_is_zero(C))
      {
         mpoly_monomial_mul_si(temp2, exp2 + 0, N, k); 
         mpoly_monomial_sub(temp2, exp, temp2, N);
         fmpz_set_mpn(t2, temp2, N);

         fmpz_divexact(temp1, C, t2);
         fmpz_add(S, S, temp1);
         fmpz_divexact(gc + gnext, temp1, poly2 + 0);

         if ((largest[1] & topbit) != 0)
         {
            x = reuse[next_free++];

            x->i = 1;
            x->j = gnext;
            x->next = NULL;

            mpoly_monomial_add(exp_list[exp_next], exp2 + N, ge + gnext*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len, N))
               exp_next--;
 
            largest[1] = gnext + 1;
         }
      }

      if (!fmpz_is_zero(S))
      {
         fmpz_set(p1 + rnext, S);
         mpoly_monomial_add(e1 + rnext*N, ge + gnext*N, exp2 + 0, N);
      } else
         rnext--;

      if (fmpz_is_zero(C))
         gnext--;

      exp_list[--exp_next] = exp;
   }

   rnext++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   fmpz_clear(t1);
   fmpz_clear(t2);
   fmpz_clear(C);
   fmpz_clear(S);
   fmpz_clear(temp1);

   flint_free(ge);
   for (i = 0; i < g_alloc; i++)
      fmpz_clear(gc + i);
   flint_free(gc);

   TMP_END;

   return rnext;
}

void fmpz_mpoly_pow_fps(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                           slong k, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2;
   ulong max = 0;
   ulong * exp2;
   int free2;

   TMP_INIT;

   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return;
   }

   if (k == 1)
   {
      fmpz_mpoly_set(poly1, poly2, ctx);

      return;
   }

   if (k == 2)
   {
      fmpz_mpoly_mul_johnson(poly1, poly2, poly2, ctx);

      return;
   }

   TMP_START;

   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);

   for (i = 0; i < ctx->n; i++)
   {
      if (max_degs2[i] > max)
         max = max_degs2[i];
   }

   if (FLINT_BIT_COUNT(max) + FLINT_BIT_COUNT(k) > sizeof(ulong)*8 || 0 > (slong) (k*max))
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_pow_fps");

   bits = FLINT_BIT_COUNT(k*max);

   exp_bits = 8;
   while (bits >= exp_bits)
      exp_bits *= 2;

   exp2 = _fmpz_mpoly_unpack_monomials(exp_bits, poly2->exps, 
                                           poly2->bits, ctx->n, poly2->length);

   free2 = exp2 != poly2->exps;

   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   if (poly2->length == 1)
   {
      fmpz_mpoly_fit_length(poly1, 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      fmpz_pow_ui(poly1->coeffs + 0, poly2->coeffs + 0, k);
      
      for (i = 0; i < N; i++)
         poly1->exps[i] = k*poly2->exps[i];

      len = 1;

      goto cleanup;
   }

   if (poly1 == poly2)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, k*(poly2->length - 1) + 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);

      if (N == 1)
      {
         len = _fmpz_mpoly_pow_fps1(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length, k);
      } else
      {
         len = _fmpz_mpoly_pow_fps(&temp->coeffs, &temp->exps, &temp->alloc,
                                     poly2->coeffs, exp2, poly2->length, k, N);
      }

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, k*(poly2->length - 1) + 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);

      if (N == 1)
      {
         len = _fmpz_mpoly_pow_fps1(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length, k);
      } else
      {
         len = _fmpz_mpoly_pow_fps(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                     poly2->coeffs, exp2, poly2->length, k, N);
      }
   }

cleanup:

   if (free2)
      flint_free(exp2);

   _fmpz_mpoly_set_length(poly1, len, ctx);

   TMP_END;
}

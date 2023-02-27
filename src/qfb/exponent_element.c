/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "qfb.h"

/*
   find which power of the base is the exponent of f
*/
ulong find_power(qfb_t f, fmpz_t n, ulong base)
{
   ulong s = 1;
      
   do
   {
      qfb_pow_ui(f, f, n, base);
      s *= base;
   } while (!qfb_is_principal_form(f, n));

   return s;
}

ulong qfb_exponent_element_stage2(qfb_t f, fmpz_t n, ulong B2_sqrt)
{
   qfb_t pow, pow2, f2;
   fmpz_t L, r;
   slong i, i2, ret = 0;
   slong depth = FLINT_BIT_COUNT(B2_sqrt) + 1;
   qfb_hash_t * qhash = qfb_hash_init(depth);
   
   fmpz_init(L);
   fmpz_init(r);
   fmpz_abs(L, n);
   fmpz_root(L, L, 4);

   qfb_init(f2);
   qfb_init(pow);
   qfb_init(pow2);
   qfb_hash_insert(qhash, f, NULL, 1, depth);

   qfb_nucomp(f2, f, f, n, L); /* large primes are odd */
   qfb_reduce(f2, f2, n);

   qfb_set(pow, f);
   
   for (i = 1; i < B2_sqrt - 1; i += 2) /* baby steps */
   {
      qfb_nucomp(pow, pow, f2, n, L);
      qfb_reduce(pow, pow, n);

      qfb_hash_insert(qhash, pow, NULL, i + 2, depth);
   }

   qfb_nucomp(pow, pow, f, n, L); /* compute f^B2_sqrt */
   qfb_reduce(pow, pow, n);

   qfb_nucomp(pow, pow, pow, n, L); /* we hash for a form or its inverse,
                                    so we can jump by f^(2x B2_sqrt) */
   qfb_reduce(pow, pow, n);
   qfb_set(pow2, pow);
   
   for(i = 2; i <= B2_sqrt; i += 2) /* giant steps */ 
   {
      i2 = qfb_hash_find(qhash, pow2, depth);
      if (i2 != -1) /* found collision */
      {
         fmpz_set_ui(r, B2_sqrt);
         fmpz_mul_ui(r, r, i);
         if (fmpz_sgn(qhash[i2].q->b) == fmpz_sgn(pow2->b))
            fmpz_sub_ui(r, r, qhash[i2].iter);
         else
            fmpz_add_ui(r, r, qhash[i2].iter);
         ret = (fmpz_size(r) > 1 ? 0 : fmpz_get_ui(r)); /* we probably should be more aggressive here */
         break;
      }

      qfb_nucomp(pow2, pow2, pow, n, L);
      qfb_reduce(pow2, pow2, n);
   }

   fmpz_clear(r);
   fmpz_clear(L);
   qfb_clear(f2);
   qfb_clear(pow);
   qfb_clear(pow2);
   qfb_hash_clear(qhash, depth);

   return ret;
}

typedef struct
{
   ulong pr;
   qfb_t pow;
   slong i;
} qfb_restart_t;

#define go_restart \
   do { \
      need_restarts = 0; \
      if (j == 0) \
      { \
         qfb_pow(f2, f2, n, exponent); \
         goto do_restart1; \
      } \
      j--; \
      qfb_set(pow, restart[j].pow); \
      qfb_pow(pow, pow, n, exponent); \
      i = restart[j].i; \
      n_primes_jump_after(iter, restart[j].pr); \
      pr = restart[j].pr; \
      goto do_restart; \
   } while (0)

int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, ulong B1, ulong B2_sqrt)
{
   slong i, j, iters = 1024, restart_inc;
   qfb_t pow, oldpow, f2;
   ulong pr, oldpr, s2, sqrt, exp;
   fmpz_t prod, L, pow2;
   int ret = 1, num_restarts = 0, need_restarts = 1, clean2 = 0;
   qfb_restart_t restart[128];
   n_primes_t iter;
   ulong hi, lo;
   double quot;
   mp_bitcnt_t bits0;

   n_primes_init(iter);

   sqrt = n_sqrt(B1);
   bits0 = FLINT_BIT_COUNT(B1);
   
   fmpz_set_ui(exponent, 1);
      
   qfb_init(f2);
   qfb_init(pow);
   qfb_init(oldpow);
   
   fmpz_init(pow2);
   fmpz_init(prod);
   fmpz_init(L);

   fmpz_abs(L, n);
   fmpz_root(L, L, 4);

   qfb_set(f2, f);

do_restart1:
   
   if (qfb_is_principal_form(f2, n))
      goto cleanup2;

   /* raise to appropriate power of 2 */
   qfb_set(oldpow, f2);
   fmpz_set_ui(pow2, 2);
   fmpz_pow_ui(pow2, pow2, bits0);
   qfb_pow(pow, oldpow, n, pow2);
   if (qfb_is_principal_form(pow, n))
   {
      ulong s = find_power(oldpow, n, 2);
      fmpz_mul_ui(exponent, exponent, s);

      goto cleanup2;
   }
   
   for (i = 0; i < 128; i++)
      qfb_init(restart[i].pow);

   clean2 = 1;

   n_prime_pi_bounds(&lo, &hi, B1);
   restart_inc = (((hi/1024 + 127)/128))*1024;
   quot = (double) B2_sqrt/(double) lo;
                      
   pr = 2;
   n_primes_jump_after(iter, 2);
   i = 0;
   
do_restart:

   /* keep raising iters by a factor of 2 until pr exceeds B1 or exponent found */
   do
   {
      /* raise to prime powers until the identity is found */
      for ( ; i < iters; )
      {
         if ((i % restart_inc) == 0 && need_restarts)
         {
            qfb_set(restart[num_restarts].pow, pow);
            restart[num_restarts].i = i;
            restart[num_restarts].pr = pr;
            num_restarts++;
         }
         
         j = FLINT_MIN(i + 1024, iters);
         qfb_set(oldpow, pow);
         oldpr = pr;
         
         for ( ; i < j; i+=4)
         {
            pr = n_primes_next(iter);
            fmpz_set_ui(prod, pr);
            if (pr < sqrt) 
            {   
               exp = bits0/FLINT_BIT_COUNT(pr);
               fmpz_pow_ui(prod, prod, exp);
            }
            pr = n_primes_next(iter);
            if (pr < sqrt) 
            {   
               exp = bits0/FLINT_BIT_COUNT(pr);
               fmpz_set_ui(pow2, pr);
               fmpz_pow_ui(pow2, pow2, exp);
               fmpz_mul(prod, prod, pow2);
            } else
               fmpz_mul_ui(prod, prod, pr);
            pr = n_primes_next(iter);
            if (pr < sqrt) 
            {   
               exp = bits0/FLINT_BIT_COUNT(pr);
               fmpz_set_ui(pow2, pr);
               fmpz_pow_ui(pow2, pow2, exp);
               fmpz_mul(prod, prod, pow2);
            } else
               fmpz_mul_ui(prod, prod, pr);
            pr = n_primes_next(iter);
            if (pr < sqrt) 
            {   
               exp = bits0/FLINT_BIT_COUNT(pr);
               fmpz_set_ui(pow2, pr);
               fmpz_pow_ui(pow2, pow2, exp);
               fmpz_mul(prod, prod, pow2);
            } else
               fmpz_mul_ui(prod, prod, pr);
            qfb_pow_with_root(pow, pow, n, prod, L);
         }

         /* identity is found, compute exponent recursively */
         if (qfb_is_principal_form(pow, n))
         {
            qfb_set(pow, oldpow);
            n_primes_jump_after(iter, oldpr);
            pr = oldpr;
            while (1)
            {
               pr = n_primes_next(iter);
         
               if (pr < sqrt)
               {
                  ulong k;
                  exp = bits0/FLINT_BIT_COUNT(pr);
                  
                  for (k = 1; k <= exp; k++)
                  {
                     qfb_pow_ui(pow, pow, n, pr);
                     if (qfb_is_principal_form(pow, n))
                     {
                        fmpz_set_ui(pow2, pr);
                        fmpz_pow_ui(pow2, pow2, k);
                        fmpz_mul(exponent, exponent, pow2);
                        for (j = 0; j < num_restarts; j++)
                        {
                           qfb_set(pow, restart[j].pow);
                           qfb_pow(pow, pow, n, exponent);
                           if (qfb_is_principal_form(pow, n))
                              break;
                        }
                        go_restart;
                     }
                  }
               } else
               {
                  qfb_pow_ui(pow, pow, n, pr);
                  if (qfb_is_principal_form(pow, n))
                  {
                     fmpz_mul_ui(exponent, exponent, pr);
                     for (j = 0; j < num_restarts; j++)
                     {
                        qfb_set(pow, restart[j].pow);
                        qfb_pow(pow, pow, n, exponent);
                        if (qfb_is_principal_form(pow, n))
                           break;
                     }
                     go_restart;
                  }
               }
            }
         }
      }
   
      /* stage 2 */
      s2 = qfb_exponent_element_stage2(pow, n, (ulong) ((double) iters * quot));
      if (s2 && n_is_prime(s2)) /* we probably should be more aggressive here */
      {
         fmpz_mul_ui(exponent, exponent, s2);
         for (j = 0; j < num_restarts; j++)
         {
            qfb_set(pow, restart[j].pow);
            qfb_pow(pow, pow, n, exponent);
            if (qfb_is_principal_form(pow, n))
               break;
         }
         go_restart;
      } 

      iters = FLINT_MIN(2*iters, hi);
   } while (pr <= B1);

   ret = 0;

cleanup2:

   if (clean2)
      for (i = 0; i < 128; i++)
         qfb_clear(restart[i].pow);

   qfb_clear(f2);
   qfb_clear(pow);
   qfb_clear(oldpow);
   fmpz_clear(pow2);
   fmpz_clear(prod);
   fmpz_clear(L);
   n_primes_clear(iter);

   return ret;
}

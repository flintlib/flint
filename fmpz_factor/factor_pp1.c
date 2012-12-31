/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_factor.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

typedef struct
{
   fmpz_t l;
   fmpz_t m;
} ppm1_struct;

typedef ppm1_struct ppm1_t[1];

void ppm1_init(ppm1_t L)
{
   fmpz_init(L->l);
   fmpz_init(L->m);
}

void ppm1_clear(ppm1_t L)
{
   fmpz_clear(L->l);
   fmpz_clear(L->m);
}

void ppm1_set(ppm1_t L, ppm1_t M)
{
   fmpz_set(L->l, M->l);
   fmpz_set(L->m, M->m);
}

void pp1_mul(ppm1_t R, ppm1_t L, ppm1_t M, const fmpz_t n, ulong c)
{
   fmpz_t t;
   fmpz_init(t);

   fmpz_mul(t, L->m, M->m);
   fmpz_mul_ui(t, t, c);
   fmpz_addmul(t, L->l, M->l);
   if (fmpz_is_odd(t))
      fmpz_add(t, t, n);
   fmpz_fdiv_q_2exp(t, t, 1);
   fmpz_mod(t, t, n);

   fmpz_mul(R->l, L->l, M->m);
   fmpz_swap(t, R->l);
   fmpz_addmul(t, L->m, M->l);
   if (fmpz_is_odd(t))
      fmpz_add(t, t, n);
   fmpz_fdiv_q_2exp(t, t, 1);
   fmpz_mod(R->m, t, n);

   fmpz_clear(t);
}

void pp1_dup(ppm1_t R, ppm1_t L, const fmpz_t n, ulong c)
{
   fmpz_t t;
   fmpz_init(t);
   
   fmpz_mul(t, L->m, L->m);
   fmpz_mul_ui(t, t, c);
   fmpz_addmul(t, L->l, L->l);
   if (fmpz_is_odd(t))
      fmpz_add(t, t, n);
   fmpz_fdiv_q_2exp(t, t, 1);
   fmpz_mod(t, t, n);

   fmpz_mul(R->l, L->l, L->m);
   fmpz_swap(t, R->l);
   fmpz_mod(R->m, t, n);

   fmpz_clear(t);
}

void pp1_pow_ui(ppm1_t R, ppm1_t L, ulong exp, const fmpz_t n, ulong c)
{
   ppm1_t P;
   ppm1_init(P);

   ppm1_set(P, L);

   while ((exp % 2) == 0)
   {
      pp1_dup(P, P, n, c);
      exp >>= 1;
   }

   ppm1_set(R, P);
   exp >>= 1;

   while (exp)
   {
      pp1_dup(P, P, n, c);
      if ((exp % 2) == 1)
         pp1_mul(R, R, P, n, c);
      exp >>= 1;
   }

   ppm1_clear(P);
}

int pp1_find_power(fmpz_t factor, ppm1_t oldL, ulong p, ulong c, const fmpz_t n)
{
   do
   {
      pp1_pow_ui(oldL, oldL, p, n, c);
      fmpz_gcd(factor, n, oldL->m);
      if (fmpz_is_zero(factor))
         return 0;
   } while (fmpz_is_one(factor));

   return 1;
}

int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, long iters, ulong c)
{
   long i, j;
   int ret = 0;
   ppm1_t L, M, oldL;
   ulong pr, oldpr;
   c = n_nextprime(c, 0);

   if (fmpz_is_even(n))
   {
      fmpz_set_ui(factor, 2);
      return 1;
   }

   iters = 128*((iters + 127)/128); /* round to a multiple of 128 iterations */

   ppm1_init(L);
   ppm1_init(oldL);
   ppm1_init(M);

   fmpz_set_ui(L->l, 1);
   fmpz_set_ui(L->m, 2);

   /* mul by various prime powers */
   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 4096, n, c); /* 2^12 */
   fmpz_gcd(factor, n, L->m);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 2, c, n);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }
   
   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 59049, n, c); /* 3^10 */
   fmpz_gcd(factor, n, L->m);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 3, c, n);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 390625, n, c); /* 5^8 */
   fmpz_gcd(factor, n, L->m);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 5, c, n);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 117649, n, c); /* 7^6 */
   fmpz_gcd(factor, n, L->m);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 7, c, n);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 14641, n, c); /* 11^4 */
   fmpz_gcd(factor, n, L->m);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 11, c, n);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   pr = 11;
   oldpr = 11;
   for (i = 0; i < iters; )
   {
      j = i + 128;
      oldpr = pr;
      for ( ; i < j; i++)
      {
         pr = n_nextprime(pr, 0);
         if (i < 131072)
            pp1_pow_ui(L, L, pr*pr, n, c);
         else
            pp1_pow_ui(L, L, pr, n, c);
      }

      fmpz_gcd(factor, n, L->m);
      if (fmpz_is_zero(factor))
         break;
      if (!fmpz_is_one(factor))
      {
         ret = 1;
         goto cleanup;
      }
   }

   if (i != iters) /* factor = 0 */
   {
      pr = oldpr;
      do
      {
         pr = n_nextprime(pr, 0);
         if (i < 131072)
            pp1_pow_ui(L, L, pr*pr, n, c);
         else
            pp1_pow_ui(L, L, pr, n, c);

         fmpz_gcd(factor, n, L->m);
         if (!fmpz_is_one(factor))
         {
            ret = (fmpz_is_zero(factor) ? 0 : 1);
            goto cleanup;
         }
      } while (1);
   }

cleanup:

   ppm1_clear(L);
   ppm1_clear(oldL);
   ppm1_clear(M);

   return ret;
}

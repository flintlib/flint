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

void ppm1_print(ppm1_t L)
{
   printf("["), fmpz_print(L->l), printf(", "), fmpz_print(L->m), printf("]");
}

void pp1_2k(ppm1_t R, ppm1_t L, const fmpz_t n, mp_limb_t ninv, fmpz_t L0, int sec)
{
   if (sec)
   {
      fmpz_mul(R->m, L->l, L->m);
      fmpz_sub(R->m, R->m, L0);
   }

   fmpz_mul(R->l, L->l, L->l);
   fmpz_sub_ui(R->l, R->l, 2);
   
   if (sec)
   {
      if (fmpz_sgn(R->m) < 0)
         fmpz_add(R->m, R->m, n);
      else
         fmpz_mod_preinv1(R->m, R->m, n, ninv);
   }

   if (fmpz_sgn(R->l) < 0)
      fmpz_add(R->l, R->l, n);
   else
      fmpz_mod_preinv1(R->l, R->l, n, ninv);
}

void pp1_2kp1(ppm1_t R, ppm1_t L, const fmpz_t n, mp_limb_t ninv, fmpz_t L0, int sec)
{
   fmpz_mul(R->l, L->l, L->m);
   fmpz_sub(R->l, R->l, L0);
   
   if (sec)
   {
      fmpz_mul(R->m, L->m, L->m);
      fmpz_sub_ui(R->m, R->m, 2);
   }

   if (fmpz_sgn(R->l) < 0)
      fmpz_add(R->l, R->l, n);
   else
      fmpz_mod_preinv1(R->l, R->l, n, ninv);

   if (sec)
   {
      if (fmpz_sgn(R->m) < 0)
         fmpz_add(R->m, R->m, n);
      else
         fmpz_mod_preinv1(R->m, R->m, n, ninv);
   }
}

void pp1_pow_ui(ppm1_t R, ppm1_t L, ulong exp, const fmpz_t n, mp_limb_t ninv)
{
   fmpz_t L0;
   
   ulong b = FLINT_BIT_COUNT(exp);
   ulong bit = (1UL<<(b - 1));

   fmpz_init(L0);
   fmpz_set(L0, L->l);

   fmpz_set(R->l, L->l);

   fmpz_mul(R->m, R->l, R->l);
   fmpz_sub_ui(R->m, R->m, 2);
   if (fmpz_sgn(R->m) < 0)
      fmpz_add(R->m, R->m, n);
   else
      fmpz_mod_preinv1(R->m, R->m, n, ninv);
   
   bit >>= 1;

   while (bit)
   {
      if (exp & bit)
         pp1_2kp1(R, R, n, ninv, L0, bit != 1);
      else
         pp1_2k(R, R, n, ninv, L0, bit != 1);

      bit >>= 1;
   }

   fmpz_clear(L0);
}

void pp1_factor(fmpz_t factor, const fmpz_t n, const fmpz_t l)
{
   fmpz_t t;

   fmpz_init(t);

   fmpz_sub_ui(t, l, 2);
   if (fmpz_sgn(t) < 0)
      fmpz_add(t, t, n);

   fmpz_gcd(factor, n, t);

   fmpz_clear(t);
}

int pp1_find_power(fmpz_t factor, ppm1_t oldL, ulong p, const fmpz_t n, mp_limb_t ninv)
{
   do
   {
      pp1_pow_ui(oldL, oldL, p, n, ninv);
      pp1_factor(factor, n, oldL->l);
      if (fmpz_is_zero(factor))
         return 0;
   } while (fmpz_is_one(factor));

   return 1;
}

int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B0, ulong c)
{
   long i, j;
   int ret = 0;
   ppm1_t L, M, oldL;
   ulong pr, oldpr, sqrt, bits0;
   mp_limb_t ninv;
   
   if (fmpz_is_even(n))
   {
      fmpz_set_ui(factor, 2);
      return 1;
   }

   sqrt = n_sqrt(B0);
   bits0 = FLINT_BIT_COUNT(B0);

   ninv = fmpz_preinv1(n);

   ppm1_init(L);
   ppm1_init(oldL);
   ppm1_init(M);

   fmpz_set_ui(L->l, c);
   fmpz_set_ui(L->m, c);
   fmpz_mul(L->m, L->m, L->m);
   fmpz_sub_ui(L->m, L->m, 2);
   if (fmpz_sgn(L->m) < 0)
      fmpz_add(L->m, L->m, n);

   /* mul by various prime powers */
   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 4096, n, ninv); /* 2^12 */
   pp1_factor(factor, n, L->l);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 2, n, ninv);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }
   
   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 59049, n, ninv); /* 3^10 */
   pp1_factor(factor, n, L->l);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 3, n, ninv);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 390625, n, ninv); /* 5^8 */
   pp1_factor(factor, n, L->l);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 5, n, ninv);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 117649, n, ninv); /* 7^6 */
   pp1_factor(factor, n, L->l);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 7, n, ninv);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   ppm1_set(oldL, L);
   pp1_pow_ui(L, L, 14641, n, ninv); /* 11^4 */
   pp1_factor(factor, n, L->l);
   if (fmpz_is_zero(factor))
   {
      ret = pp1_find_power(factor, oldL, 11, n, ninv);
      goto cleanup;
   }
   if (!fmpz_is_one(factor))
   {
      ret = 1;
      goto cleanup;
   }

   pr = 11;
   oldpr = 11;
   for (i = 0; pr < B0; )
   {
      j = i + 1024;
      oldpr = pr;
      for ( ; i < j; i++)
      {
         pr = n_nextprime(pr, 0);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            pp1_pow_ui(L, L, n_pow(pr, exp), n, ninv);
         } else
            pp1_pow_ui(L, L, pr, n, ninv);
      }
      
      pp1_factor(factor, n, L->l);
      if (fmpz_is_zero(factor))
         break;
      if (!fmpz_is_one(factor))
      {
         ret = 1;
         goto cleanup;
      }
   }

   if (pr < B0) /* factor = 0 */
   {
      pr = oldpr;
      do
      {
         pr = n_nextprime(pr, 0);
         if (pr < sqrt)
         {
            ulong bits = FLINT_BIT_COUNT(pr);
            ulong exp = bits0 / bits;
            pp1_pow_ui(L, L, n_pow(pr, exp), n, ninv);
         } else
            pp1_pow_ui(L, L, pr, n, ninv);

         pp1_factor(factor, n, L->l);
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

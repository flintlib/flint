/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void _fmpz_nm1_trial_factors(const fmpz_t n, mp_ptr pm1, slong * num_pm1, ulong limit)
{
   slong i, num;
   ulong ppi, p;
   const ulong * primes;
   const double * pinv;
   
   *num_pm1 = 0;

   /* number of primes multiplied that will fit in a word */
  
   num = FLINT_BITS/FLINT_BIT_COUNT(limit);

   /* compute remainders of n mod p for primes p up to limit (approx.) */

   n_prime_pi_bounds(&ppi, &ppi, limit); /* precompute primes */
   primes = n_primes_arr_readonly(ppi + FLINT_BITS);
   pinv = n_prime_inverses_arr_readonly(ppi + FLINT_BITS);
   
   while (primes[0] < limit)
   {
      /* multiply batch of primes */
      
      p = primes[0];
      for (i = 1; i < num; i++)
         p *= primes[i];

      /* multi-modular reduction */

      p = fmpz_tdiv_ui(n, p);

      /* check for factors */
      for (i = 0; i < num; i++)
      {
         ulong r = n_mod2_precomp(p, primes[i], pinv[i]);

         if (r == 1) /* n - 1 = 0 mod p */
            pm1[(*num_pm1)++] = primes[i];
      }

      /* get next batch of primes */
      primes += num;
      pinv += num;
   }
}

int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pm1, slong num_pm1)
{
   slong i, d, bits;
   ulong a;
   fmpz_t g, q, r, pow, pow2, ex, c, p;
   fmpz_factor_t fac;
   int res = 0, fac_found;
      
   fmpz_init(p);
   fmpz_init(q);
   fmpz_init(r);
   fmpz_init(g);
   fmpz_init(pow);
   fmpz_init(pow2);
   fmpz_init(c);
   fmpz_init(ex);
   fmpz_factor_init(fac);

   fmpz_sub_ui(R, n, 1); /* start with n - 1 */

   bits = fmpz_bits(R);

   for (i = 0; i < num_pm1; i++)
   {
      fmpz_set_ui(p, pm1[i]);
      d = fmpz_remove(R, R, p);
      _fmpz_factor_append_ui(fac, pm1[i], d);
   }

   srand(time(NULL));

   if (!fmpz_is_probabprime_BPSW(R))  
   {
      if (bits > 150 && (fac_found = fmpz_factor_pp1(p, R, bits + 1000, bits/20 + 1000, rand()%100 + 3)
                    && fmpz_is_prime(p)))
      {
         d = fmpz_remove(R, R, p);
         _fmpz_factor_append(fac, p, d);

         if (fmpz_is_probabprime_BPSW(R)) /* fast test first */
         {
            if (fmpz_is_prime(R) == 1)
            {
               _fmpz_factor_append(fac, R, 1);
               fmpz_set_ui(R, 1);
            }
         } 
      }
   } else
   {
      if (fmpz_is_prime(R) == 1)
      {
         _fmpz_factor_append(fac, R, 1);
         fmpz_set_ui(R, 1);
      }
   } 

   /* compute product F of found primes */
   fmpz_set_ui(F, 1);
   for (i = 0; i < fac->num; i++)
   {
      if (fac->exp[i] == 1)
         fmpz_mul(F, F, fac->p + i);
      else
      {
         fmpz_pow_ui(pow, fac->p + i, fac->exp[i]);
         fmpz_mul(F, F, pow);
      }
   }
   
   for (a = 2; ; a++)
   {
      /* compute a^((n-1)/F) mod n */
      fmpz_set_ui(pow, a);
      fmpz_powm(pow, pow, R, n);
      
      /* check a^(n-1) = 1 mod n */
      fmpz_powm(pow2, pow, F, n);
      if (!fmpz_is_one(pow2))
      {
         res = 0;
         goto cleanup;
      }

      fmpz_set_ui(c, 1);

      /* find values a^((n-1)/q) - 1 for each prime q dividing F */
      for (i = 0; i < fac->num; i++)
      {
         fmpz_tdiv_q(ex, F, fac->p + i);
         fmpz_powm(pow2, pow, ex, n);
         fmpz_sub_ui(pow2, pow2, 1);
         if (fmpz_sgn(pow2) < 0)
            fmpz_add(pow2, pow2, n);

         if (!fmpz_is_zero(pow2))
         {
            fmpz_mul(c, c, pow2);
            fmpz_mod(c, c, n);
         } else
            break;
      }

      if (i == fac->num) /* found valid base a */
         break;
   }

   /* check for factors of n */
   fmpz_gcd(g, n, c);
   res = fmpz_is_one(g);
   
cleanup:

   fmpz_factor_clear(fac);
   fmpz_clear(pow);
   fmpz_clear(pow2);
   fmpz_clear(c);
   fmpz_clear(ex);
   fmpz_clear(p);
   fmpz_clear(q);
   fmpz_clear(r);
   fmpz_clear(g);
   
   return res;
}

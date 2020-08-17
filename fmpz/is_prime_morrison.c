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

void _fmpz_np1_trial_factors(const fmpz_t n, mp_ptr pp1, slong * num_pp1, ulong limit)
{
   slong i, num;
   ulong ppi, p;
   const ulong * primes;
   const double * pinv;
   
   *num_pp1 = 0;

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

         if (r == primes[i] - 1) /* n + 1 = 0 mod p */
            pp1[(*num_pp1)++] = primes[i];
      }

      /* get next batch of primes */
      primes += num;
      pinv += num;
   }
}

int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, mp_ptr pp1, slong num_pp1)
{
   slong i, d, bits;
   mp_limb_t a, b;
   fmpz_t g, q, r, ex, c, D, Dinv, A, B, Ukm, Ukm1, Um, Um1, Vm, Vm1, p;
   fmpz_factor_t fac;
   int res = 0, fac_found;

   fmpz_init(D);
   fmpz_init(Dinv);
   fmpz_init(A);
   fmpz_init(B);
   fmpz_init(p);
   fmpz_init(q);
   fmpz_init(r);
   fmpz_init(g);
   fmpz_init(c);
   fmpz_init(ex);
   fmpz_init(Um);
   fmpz_init(Um1);
   fmpz_init(Ukm);
   fmpz_init(Ukm1);
   fmpz_init(Vm);
   fmpz_init(Vm1);
   fmpz_factor_init(fac);

   fmpz_add_ui(R, n, 1); /* start with n + 1 */
   
   bits = fmpz_bits(R);

   for (i = 0; i < num_pp1; i++)
   {
      fmpz_set_ui(p, pp1[i]);
      d = fmpz_remove(R, R, p);
      _fmpz_factor_append_ui(fac, pp1[i], d);
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
         fmpz_pow_ui(ex, fac->p + i, fac->exp[i]);
         fmpz_mul(F, F, ex);
      }
   }
   
   /* Want D = A^2 - 4B where A = a, B = b, such that (D/n) = -1 */
   for (b = 1; ; b++)
   {
      fmpz_set_ui(B, b);
      fmpz_gcd(g, B, n);

      if (fmpz_equal(g, n)) /* need gcd(n, b) = 1 */
         continue;

      if (!fmpz_is_one(g)) /* found a factor of n */
      {
         res = 0;
         goto cleanup;
      }

      a = 2;
      do {
         a++;
         fmpz_set_ui(A, a);
         fmpz_mul_ui(D, A, a);
         fmpz_sub_ui(D, D, 4*b);
      } while (fmpz_jacobi(D, n) != -1);

      fmpz_invmod(Dinv, D, n);

      /* compute U((n+1)/F) mod n */
      fmpz_lucas_chain_full(Vm, Vm1, A, B, R, n);
      fmpz_lucas_chain_VtoU(Um, Um1, Vm, Vm1, A, B, Dinv, n);

      /* check U(n+1) = 0 mod n */
      fmpz_lucas_chain_mul(Ukm, Ukm1, Um, Um1, A, B, F, n);
      if (!fmpz_is_zero(Ukm))
      {
         res = 0;
         goto cleanup;
      }

      fmpz_set_ui(c, 1);

      /* find values U((n+1)/q) for each prime q dividing F */
      for (i = 0; i < fac->num; i++)
      {
         fmpz_tdiv_q(ex, F, fac->p + i);
         fmpz_lucas_chain_mul(Ukm, Ukm1, Um, Um1, A, B, ex, n);

         if (!fmpz_is_zero(Ukm))
         {
            fmpz_mul(c, c, Ukm);
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
   fmpz_clear(D);
   fmpz_clear(Dinv);
   fmpz_clear(A);
   fmpz_clear(B);
   fmpz_clear(c);
   fmpz_clear(ex);
   fmpz_clear(p);
   fmpz_clear(q);
   fmpz_clear(r);
   fmpz_clear(g);
   fmpz_clear(Um);
   fmpz_clear(Um1);
   fmpz_clear(Ukm);
   fmpz_clear(Ukm1);
   fmpz_clear(Vm);
   fmpz_clear(Vm1);

   return res;
}

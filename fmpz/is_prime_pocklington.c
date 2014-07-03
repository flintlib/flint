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

    Copyright (C) 2014 William Hart
   
******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_is_prime_pocklington(fmpz_t F, const fmpz_t n, ulong limit)
{
   n_primes_t iter;
   slong num, i, count, num2;
   mp_limb_t p = 0, a;
   fmpz * vec, * vec2;
   fmpz_t nm1; /* n - 1 */
   fmpz_t g, q, r, pow, pow2, ex, c;
   fmpz_factor_t fac;
   int res = 0;

   /* number of primes multiplied that's about the same size as n */
   num = fmpz_bits(n)/FLINT_BIT_COUNT(limit) + 1;
   /* round down to power of 2 */
   num = WORD(1) << (FLINT_BIT_COUNT(num) - 1);

   vec = _fmpz_vec_init(num + (num + 1)/2);
   vec2 = vec + num;

   n_primes_init(iter);

   fmpz_init(q);
   fmpz_init(r);
   fmpz_init(g);
   fmpz_init(pow);
   fmpz_init(pow2);
   fmpz_init(c);
   fmpz_init(ex);
   fmpz_init(nm1);
   fmpz_sub_ui(nm1, n, 1);
   fmpz_factor_init(fac);

   while (p < limit)
   {
      /* get next batch of primes */
      for (count = 0; count < num && p < limit; count++)
      {
         p = n_primes_next(iter);
         fmpz_set_ui(vec + count, p);
      }
 
      /* multiply batch */
      num2 = count;
      if (num2 > 1) /* vec -> vec2, preserve vec */
      {
         for (i = 0; i < num2/2; i++)
            fmpz_mul(vec2 + i, vec + 2*i, vec + 2*i + 1);
         
         if (num2 & 1)
            fmpz_set(vec2 + i, vec + 2*i);

         num2 = (num2 + 1)/2;
      }

      while (num2 > 1) /* vec2 -> vec2 */
      {
         for (i = 0; i < num2/2; i++)
            fmpz_mul(vec2 + i, vec2 + 2*i, vec2 + 2*i + 1);
         
         if (num2 & 1)
            fmpz_set(vec2 + i, vec2 + 2*i);

         num2 = (num2 + 1)/2;
      }

      /* gcd product with n - 1 */
      fmpz_gcd(g, nm1, vec2 + 0);

      if (!fmpz_is_one(g)) /* we have some factors */
      {
         for (i = 0; i < count; i++)
         {
            if (fmpz_equal(g, vec + i))
            {
               _fmpz_factor_append_ui(fac, fmpz_get_ui(vec + i), 1);
               break;
            }

            fmpz_tdiv_qr(q, r, g, vec + i);

            if (fmpz_is_zero(r))
            {
               fmpz_set(g, q);
               _fmpz_factor_append_ui(fac, fmpz_get_ui(vec + i), 1);
            }
         }
      }
   }

   /* compute product F of found primes */
   fmpz_set_ui(F, 1);
   for (i = 0; i < fac->num; i++)
      fmpz_mul(F, F, fac->p + i);
   
   for (a = 2; ; a++)
   {
      /* compute a^((n-1)/F) mod n */
      fmpz_tdiv_q(ex, nm1, F);
      fmpz_set_ui(pow, a);
      fmpz_powm(pow, pow, ex, n);
      
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
   fmpz_clear(q);
   fmpz_clear(r);
   fmpz_clear(g);
   fmpz_clear(nm1);
   _fmpz_vec_clear(vec, num + (num + 1)/2);
   n_primes_clear(iter);

   return res;
}

/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_factor_trial_gcd(n_factor_t * factors, mp_limb_t n)
{
   unsigned int exp;
   mp_limb_t p, n1, b, c;
   double ppre;
   ulong i;

   n_compute_primes(FLINT_NUM_FACTOR_GCD_PRIMES/3);

   n1 = n_factor_trial(factors, n, 200);
   if (n1 == 1UL) return 1UL;

   mp_limb_t a = mpn_mod_1(prime_prod, prime_prod_n, n1);

   if (a != 0UL)
   {
      b = n_gcd(n1, a);
      if (b == 1UL) return n1;
    
      c = n1/b;
      if ((b < FLINT_FACTOR_GCD_LIMIT) && (n_is_prime(b)))
      {
         exp = n_remove(&n1, b);
         n_factor_insert(factors, b, exp);
         return c;
      } else if ((c < FLINT_FACTOR_GCD_LIMIT) && (n_is_prime(c)))
      {
         exp = n_remove(&n1, c);
         n_factor_insert(factors, c, exp);          
      }
   }
 
   return n_factor_trial_range(factors, n1, 200, FLINT_NUM_FACTOR_GCD_PRIMES/3);      
}

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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_factor_trial_partial(n_factor_t * factors, mp_limb_t n, mp_limb_t * prod, ulong num_primes, mp_limb_t limit)
{
   unsigned int exp;
   mp_limb_t p;
   double ppre;
   ulong i;

   (*prod) = 1;
   n_compute_primes(num_primes);

   for (i = 0; i < num_primes; i++)
   {
      p = flint_primes[i];
      if (p*p > n) break;
      ppre = flint_prime_inverses[i];
      exp = n_remove2_precomp(&n, p, ppre);
      if (exp) 
      {
         n_factor_insert(factors, p, exp);
         (*prod) *= n_pow(p, exp);
         if (*prod > limit) break;
      }
   }
       
   return n;
}

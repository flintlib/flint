/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_factor_trial_partial(n_factor_t * factors, mp_limb_t n, mp_limb_t * prod, ulong num_primes, mp_limb_t limit)
{
   unsigned int exp;
   mp_limb_t p;
   double ppre;
   ulong i;
   const mp_limb_t * primes;
   const double * inverses;

   (*prod) = 1;
   primes = n_primes_arr_readonly(num_primes);
   inverses = n_prime_inverses_arr_readonly(num_primes);

   for (i = 0; i < num_primes; i++)
   {
      p = primes[i];
      if (p*p > n) break;
      ppre = inverses[i];
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

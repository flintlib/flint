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

mp_limb_t n_factor_trial_range(n_factor_t * factors, mp_limb_t n, ulong start, ulong num_primes)
{
   unsigned int exp;
   mp_limb_t p;
   double ppre;
   ulong i;
   const mp_limb_t * primes;
   const double * inverses;

   primes = n_primes_arr_readonly(num_primes);
   inverses = n_prime_inverses_arr_readonly(num_primes);

   for (i = start; i < num_primes; i++)
   {
      p = primes[i];
      if (p*p > n) break;
      ppre = inverses[i];
      exp = n_remove2_precomp(&n, p, ppre);
      if (exp) n_factor_insert(factors, p, exp);
   }
       
   return n;
}

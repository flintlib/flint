/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_factor_trial(n_factor_t * factors, ulong n, ulong num_primes)
{
    return n_factor_trial_range(factors, n, UWORD(0), num_primes);
}

ulong n_factor_trial_partial(n_factor_t * factors, ulong n, ulong * prod, ulong num_primes, ulong limit)
{
   unsigned int exp;
   ulong p;
   double ppre;
   ulong i;
   const ulong * primes;
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

ulong n_factor_trial_range(n_factor_t * factors, ulong n, ulong start, ulong num_primes)
{
   unsigned int exp;
   ulong p;
   double ppre;
   ulong i;
   const ulong * primes;
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

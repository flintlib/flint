/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

int is_prime2(mp_limb_t n, int proved)
{
   if (proved) return n_is_prime(n);
   else return n_is_probabprime(n);
}

mp_limb_t n_factor_partial(n_factor_t * factors, mp_limb_t n, mp_limb_t limit, int proved)
{
   ulong factor_arr[FLINT_MAX_FACTORS_IN_LIMB];
   ulong exp_arr[FLINT_MAX_FACTORS_IN_LIMB];
   ulong factors_left;
   ulong exp;
   mp_limb_t cofactor, factor, cutoff, prod;

   cofactor = n_factor_trial_partial(factors, n, &prod, FLINT_FACTOR_TRIAL_PRIMES, limit);
   if (prod > limit) return cofactor;

   if ( cofactor == 1 )
   {
      return cofactor;
   }

   if (is_prime2(cofactor, proved))
   {
      n_factor_insert(factors, cofactor, UWORD(1));
      return 1;
   }

   factor_arr[0] = cofactor;
   factors_left = 1;
   exp_arr[0] = 1;

   cutoff = FLINT_FACTOR_TRIAL_CUTOFF;

   while (factors_left > 0 && prod <= limit)
   {
      factor = factor_arr[factors_left - 1];

      if (factor >= cutoff)
		{
	      if ((cofactor = n_factor_power235(&exp, factor)))
         {
            exp_arr[factors_left - 1] *= exp;
            factor_arr[factors_left - 1] = factor = cofactor;
         }

         if ((factor >= cutoff) && !is_prime2(factor, proved))
		   {
		      if ((
#if FLINT64
                 (factor < FLINT_FACTOR_ONE_LINE_MAX) &&
#endif
                 (cofactor = n_factor_one_line(factor, FLINT_FACTOR_ONE_LINE_ITERS)))
              || (cofactor = n_factor_SQUFOF(factor, FLINT_FACTOR_SQUFOF_ITERS)))
				{
					exp_arr[factors_left] = exp_arr[factors_left - 1];
               factor_arr[factors_left] = cofactor;
               factor_arr[factors_left - 1] /= cofactor;
               factors_left++;
				} else
				{
               flint_throw(FLINT_ERROR, "Error (n_factor_partial). Failed to factor %wd.\n", n);
				}
         } else
			{
				n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
            prod *= n_pow(factor, exp_arr[factors_left - 1]);
            factors_left--;
			}
		} else
		{
			n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
         prod *= n_pow(factor, exp_arr[factors_left - 1]);
         factors_left--;
		}
   }

   return n/prod;
}

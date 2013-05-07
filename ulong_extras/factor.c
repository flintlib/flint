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

#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int is_prime(mp_limb_t n, int proved)
{
    return proved ? n_is_prime(n) : n_is_probabprime(n);
}

void n_factor(n_factor_t * factors, mp_limb_t n, int proved)
{
   ulong factor_arr[FLINT_MAX_FACTORS_IN_LIMB];
   ulong exp_arr[FLINT_MAX_FACTORS_IN_LIMB];
   ulong factors_left;
   ulong exp;
   mp_limb_t cofactor, factor, cutoff;

   cofactor = n_factor_trial(factors, n, FLINT_FACTOR_TRIAL_PRIMES);
   if (cofactor == 1UL) return;
   if (is_prime(cofactor, proved)) 
   {
      n_factor_insert(factors, cofactor, 1UL);
      return;
   }

   factor_arr[0] = cofactor;
   factors_left = 1;
   exp_arr[0] = 1;

   cutoff = flint_primes[FLINT_FACTOR_TRIAL_PRIMES - 1]*flint_primes[FLINT_FACTOR_TRIAL_PRIMES - 1];

   while (factors_left > 0)
   {
      factor = factor_arr[factors_left - 1];

      if (factor >= cutoff)
		{
	      if ((cofactor = n_factor_power235(&exp, factor)))
         {
            exp_arr[factors_left - 1] *= exp;
            factor_arr[factors_left - 1] = factor = cofactor;
         }
           
         if ((factor >= cutoff) && !is_prime(factor, proved))
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
               printf("Exception (n_factor). Failed to factor %ld.\n", n);
               abort();
				}
         } else
			{
				n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
            factors_left--;
			}
		} else
		{
			n_factor_insert(factors, factor, exp_arr[factors_left - 1]);
         factors_left--;
		}
   } 
}

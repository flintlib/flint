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

mp_limb_t n_factor_partial1(n_factor_t * factors, mp_limb_t n)
{
   mp_limb_t n1, n2;
   mp_limb_t cofactor = n_factor_trial_gcd(factors, n);

   if (cofactor == 1UL) return 1UL;

   if (!n_is_prime(cofactor))
   {
      n1 = n_factor_one_line(cofactor, FLINT_FACTOR_PARTIAL1_CUTOFF);
      if (!n1) return cofactor;

      n2 = cofactor/n1;

      if (!n_is_prime(n2)) 
      {
         if (!n_is_prime(n1)) return cofactor;

         factors->p[factors->num] = n1;
         factors->exp[factors->num++] = n_remove(&cofactor, n1);

         return cofactor;
      } else
      {
         factors->p[factors->num] = n2;
         factors->exp[factors->num++] = n_remove(&cofactor, n2);

         if (!n_is_prime(cofactor)) return cofactor;

         factors->p[factors->num] = cofactor;
         factors->exp[factors->num++] = 1;

         return 1UL;
      }
   } else
   {
      factors->p[factors->num] = cofactor;
      factors->exp[factors->num++] = 1;

      return 1UL;
   }
}

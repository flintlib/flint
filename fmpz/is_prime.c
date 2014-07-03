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
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_is_prime(const fmpz_t n)
{
   int res = 0;
   
   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;
   
   if (fmpz_is_even(n))
      return (fmpz_cmp_ui(n, 2) == 0);

   {
      double logd = log(fmpz_get_d(n));
      ulong limit = (ulong) (logd*logd*logd/10.0) + 2;

      fmpz_t F, F2;

      fmpz_init(F);
      fmpz_init(F2);

      res = fmpz_is_prime_pocklington(F, n, limit);

      if (res == 1)
      {
         fmpz_mul(F2, F, F);
         if (fmpz_cmp(F2, n) < 0)
            res = -1;
      }

      fmpz_clear(F);
      fmpz_clear(F2);
   }

   return res;
      
}

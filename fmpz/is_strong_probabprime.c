/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t base)
{
   fmpz_t a, nm1, t, y;
   int res = 0;

   if (fmpz_cmp_ui(n, 1) <= 0)
      return 0;

   fmpz_init(a);
   fmpz_init(t);
   fmpz_init(nm1);
   
   fmpz_sub_ui(nm1, n, 1);

   if (fmpz_cmp(base, n) >= 0)
      fmpz_mod(a, base, n);
   else
      fmpz_set(a, base);

   if (fmpz_is_one(a) || fmpz_equal(a, nm1) || fmpz_is_zero(a))
      res = 1;
   else
   {
      slong s = 0;

      fmpz_init(y);
      s = fmpz_val2(nm1);

      fmpz_tdiv_q_2exp(t, nm1, s);
      fmpz_powm(y, a, t, n);

      if (fmpz_is_one(y))
         res = 1;
      else
      {
         for (s--; s > 0 && !fmpz_equal(y, nm1); s--) 
         {
            fmpz_mul(t, y, y);
            fmpz_mod(y, t, n);
         }

         res = fmpz_equal(y, nm1);
      }

      fmpz_clear(y);
   }

   fmpz_clear(nm1);
   fmpz_clear(a);
   fmpz_clear(t);

   return res;
}

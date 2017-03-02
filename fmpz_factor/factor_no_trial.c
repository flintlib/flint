/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "qsieve.h"

void
fmpz_factor_no_trial(fmpz_factor_t factor, const fmpz_t n)
{
   int exp, i;

   if (fmpz_is_prime(n))
      _fmpz_factor_append(factor, n, 1);
   else
   {
      fmpz_t root;

      fmpz_init(root);

      exp = fmpz_is_perfect_power(root, n);

      if (exp != 0)
      {
         fmpz_factor_t fac;

         fmpz_factor_init(fac);

         fmpz_factor_no_trial(fac, root);

         _fmpz_factor_concat(factor, fac, exp);

         fmpz_factor_clear(fac);
      } else
      {
         fmpz_factor_t fac, fac2;

         fmpz_factor_init(fac);

         /* insert call to ecm here */

         qsieve_factor(fac, n);

         for (i = 0; i < fac->num; i++)
         {
            fmpz_factor_init(fac2);

            fmpz_factor_no_trial(fac2, fac->p + i);

            _fmpz_factor_concat(factor, fac2, fac->exp[i]);
 
            fmpz_factor_clear(fac2);
         }

         fmpz_factor_clear(fac);
      }
   }
}

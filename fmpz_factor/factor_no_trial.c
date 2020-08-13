/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "thread_support.h"

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
         fmpz_factor_t fac, fac2, fac3;
	 slong bits = fmpz_sizeinbase(n, 2);
         int done;

         fmpz_factor_init(fac3);

	 done = fmpz_factor_smooth(fac3, n, FLINT_MAX(bits/3 - 17, 2), 1);

         if (!done)
	 {
            fmpz_t n2;
	    slong exp2;

	    fmpz_init(n2);
	    
	    /* take out cofactor and factor it */
	    fmpz_set(n2, fac3->p + fac3->num - 1);
	    exp = fac3->exp[fac3->num - 1];
	    fac3->exp[fac3->num - 1] = 0;
	    fac3->num--;

	    fmpz_factor_init(fac);

	    /* qsieve can't factor perfect powers */
	    exp2 = fmpz_is_perfect_power(root, n2);

	    if (exp2)
                _fmpz_factor_append(fac, root, exp2);
	    else
	        qsieve_factor(fac, n2);

            for (i = 0; i < fac->num; i++)
            {
	       fmpz_factor_init(fac2);

               fmpz_factor_no_trial(fac2, fac->p + i);

               _fmpz_factor_concat(fac3, fac2, exp*fac->exp[i]);
 
               fmpz_factor_clear(fac2);
            }

            fmpz_factor_clear(fac);

	    fmpz_clear(n2);
	 }

	 _fmpz_factor_concat(factor, fac3, 1);

	 fmpz_factor_clear(fac3);
      }

      fmpz_clear(root);
   }
}

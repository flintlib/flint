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
#include "fmpz_mpoly.h"

void mpoly_max_degrees1(ulong * max_degs, const ulong * exps,
                                                slong len, slong bits, slong n)
{
   slong i, j, k = FLINT_BITS/bits, shift;
   ulong mask;

   shift = (k - n)*bits;
   mask = bits == FLINT_BITS ? ~UWORD(0) : (UWORD(1) << bits) - UWORD(1);

   for (i = 0; i < len; i++) /* for each exponent vector */
   {
      ulong v = (exps[i] >> shift); /* shift unused exponents */

      for (j = 0; j < n; j++) /* for each exponent */
      {
         ulong ex = (v & mask);

         if (ex > max_degs[j])
            max_degs[j] = ex;

         v >>= bits;
      }      
   }
}

void mpoly_max_degrees(ulong * max_degs, const ulong * poly_exps,
                                                slong len, slong bits, slong n)
{
   slong i, j;
   ulong * exps;
   slong N;
   TMP_INIT;
   
   for (i = 0; i < n; i++)
      max_degs[i] = 0;

   N = (bits*n - 1)/FLINT_BITS + 1;
   
   if (N == 1)
   {
       mpoly_max_degrees1(max_degs, poly_exps, len, bits, n);

       return;
   }

   TMP_START;
     
   exps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   
   for (i = 0; i < len; i++)
   {
      mpoly_get_monomial(exps, poly_exps + i*N, bits, n, 0, 0);

      for (j = 0; j < n; j++)
      {
         if (exps[n - j - 1] > max_degs[j])
            max_degs[j] = exps[n - j - 1];
      }
   }

   TMP_END;
}

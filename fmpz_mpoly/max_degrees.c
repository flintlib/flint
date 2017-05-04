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

void _fmpz_mpoly_max_degrees1(ulong * max_degs, const ulong * exps,
                              slong len, slong bits, slong n, int deg, int rev)
{
   slong i, j, k = FLINT_BITS/bits, shift;
   ulong mask;

   for (i = 0; i < n; i++)
      max_degs[i] = 0;

   shift = (k - n)*bits;
   mask = bits == FLINT_BITS ? ~UWORD(0) : (UWORD(1) << bits) - UWORD(1);

   for (i = 0; i < len; i++) /* for each exponent vector */
   {
      ulong v = (exps[i] >> shift); /* shift unused exponents */

      if (rev)
      {
         for (j = 0; j < n - deg; j++) /* for each exponent except degree */
         {
            ulong ex = (v & mask);

            if (ex > max_degs[j + deg])
               max_degs[j + deg] = ex;

            v >>= bits;
         }
      } else
      {
         for (j = 0; j < n - deg; j++) /* for each exponent except degree */
         {
            ulong ex = (v & mask);

            if (ex > max_degs[n - j - 1])
               max_degs[n - j - 1] = ex;

            v >>= bits;
         }
      }      

      if (deg)
      {
         if (v > max_degs[0])
            max_degs[0] = v;
      }
   }
}

void _fmpz_mpoly_max_degrees(ulong * max_degs, const ulong * poly_exps,
                              slong len, slong bits, slong n, int deg, int rev, slong N)
{
   slong i, j;
   ulong * exps;
   TMP_INIT;
   
   if (N == 1)
   {
       _fmpz_mpoly_max_degrees1(max_degs, poly_exps, len, bits, n, deg, rev);

       return;
   }

   TMP_START;
     
   exps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   
   for (i = 0; i < len; i++)
   {
      ulong s = 0;

      mpoly_get_monomial(exps, poly_exps + i*N, bits, n, deg, rev);

      for (j = 0; j < n - deg; j++)
      {
         s += exps[j];

         if (exps[j] > max_degs[j + deg])
            max_degs[j + deg] = exps[j];
      }

      if (deg && s > max_degs[0])
         max_degs[0] = s;
   }

   TMP_END;
}

void fmpz_mpoly_max_degrees(ulong * max_degs, const fmpz_mpoly_t poly,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;
   
   slong N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   degrev_from_ord(deg, rev, ctx->ord);

   _fmpz_mpoly_max_degrees(max_degs, poly->exps, poly->length,
                                              poly->bits, ctx->n, deg, rev, N);
}

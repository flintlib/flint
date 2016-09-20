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

ulong * _fmpz_mpoly_max_degrees1(const ulong * exps, slong len, slong bits, slong n, int deg, int rev)
{
   slong i, j, k = FLINT_BITS/bits, shift;
   ulong * maxdegs = (ulong *) flint_malloc(n*sizeof(ulong));
   ulong mask;

   for (i = 0; i < n; i++)
      maxdegs[i] = 0;

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

            if (ex > maxdegs[j + deg])
               maxdegs[j + deg] = ex;

            v >>= bits;
         }
      } else
      {
         for (j = 0; j < n - deg; j++) /* for each exponent except degree */
         {
            ulong ex = (v & mask);

            if (ex > maxdegs[n - j - 1])
               maxdegs[n - j - 1] = ex;

            v >>= bits;
         }
      }      

      if (deg)
      {
         if (v > maxdegs[0])
            maxdegs[0] = v;
      }
   }

   return maxdegs;
}

ulong * _fmpz_mpoly_max_degrees(const ulong * poly_exps,
                              slong len, slong bits, slong n, int deg, int rev)
{
   slong i, j, N;
   ulong * exps, * max_exps;
   TMP_INIT;
   
   TMP_START;
     
   N = (bits*n - 1)/FLINT_BITS + 1;

   exps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   max_exps = (ulong *) flint_calloc(n, sizeof(ulong));

   for (i = 0; i < len; i++)
   {
      _fmpz_mpoly_get_monomial(exps, poly_exps + i*N, bits, n, deg, rev);

      for (j = 0; j < n; j++)
      {
         if (exps[j] > max_exps[j])
            max_exps[j] = exps[j];
      }
   }

   TMP_END;

   return max_exps;
}

ulong * fmpz_mpoly_max_degrees(const fmpz_mpoly_t poly,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;
   
   slong N = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   degrev_from_ord(deg, rev, ctx->ord);

   if (N == 1)
      return _fmpz_mpoly_max_degrees1(poly->exps, poly->length,
                                                  poly->bits, ctx->n, deg, rev);
   else
      return _fmpz_mpoly_max_degrees(poly->exps, poly->length,
                                                  poly->bits, ctx->n, deg, rev);
   
   return NULL;
}

/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;
   slong i;

   if (poly->length != 1)
      return 0;

   if (!fmpz_is_one(poly->coeffs + 0))
      return 0;

   degrev_from_ord(deg, rev, ctx->ord);

   if (deg)
      return (poly->exps[0] >> (FLINT_BITS - poly->bits)) == 1;
   else
   {
      slong m = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

      for (i = 0; i < m; i++)
      {
         if (poly->exps[i] != 0)
         {
            ulong s = 1;

            while (s != 0 && poly->exps[i] != s)
               s <<= poly->bits;

            if (poly->exps[i] != s)
               return 0;

            for (i++; i < m; i++)
            {
               if (poly->exps[i] != 0)
                  return 0;
            }

            return 1;
         }
      }

      return 0;
   }
}

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

void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                     const fmpz * poly2, const ulong * exps2, slong n, slong m)
{
   slong i;

   if (poly1 != poly2)
   {
      for (i = 0; i < n; i++)
         fmpz_set(poly1 + i, poly2 + i);
   }

   if (exps1 != exps2)
   {
      for (i = 0; i < n*m; i++)
         exps1[i] = exps2[i];
   }
}

void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   slong m = (poly2->bits*ctx->n - 1)/FLINT_BITS + 1;

   fmpz_mpoly_fit_length(poly1, poly2->length, ctx);
   fmpz_mpoly_fit_bits(poly1, poly2->bits, ctx);

   _fmpz_mpoly_set(poly1->coeffs, poly1->exps,
                   poly2->coeffs, poly2->exps, poly2->length, m);

   _fmpz_mpoly_set_length(poly1, poly2->length, ctx);
   poly1->bits = poly2->bits;
}

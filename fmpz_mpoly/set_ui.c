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

void fmpz_mpoly_set_ui(fmpz_mpoly_t poly, ulong c, const fmpz_mpoly_ctx_t ctx)
{
   slong m, i;

   if (c == 0)
   {
      _fmpz_mpoly_set_length(poly, 0, ctx);
      return;
   }

   fmpz_mpoly_fit_length(poly, 1, ctx);

   fmpz_set_ui(poly->coeffs, c);

   m = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;

   for (i = 0; i < m; i++)
      poly->exps[i] = 0;

   _fmpz_mpoly_set_length(poly, 1, ctx);
}

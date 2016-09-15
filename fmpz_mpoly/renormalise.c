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

void _fmpz_mpoly_renormalise(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   slong i = 0, j;
   slong m;
   fmpz * tmp;

   TMP_INIT;

   while (i < poly->length && fmpz_is_zero(poly->coeffs + i))
      i++;

   if (i != 0)
   {
      m = (poly->bits*ctx->n - 1)/FLINT_BITS + 1;;

      TMP_START;

      tmp = (fmpz *) TMP_ALLOC(i*sizeof(fmpz));

      for (j = 0; j < i; j++)
         tmp[j] = poly->coeffs[j];

      for (j = i; j < poly->length; j++)
         poly->coeffs[j - i] = poly->coeffs[j];

      for (j = 0; j < i; j++)
         poly->coeffs[j + i] = tmp[j];

      TMP_END;

      for (j = i*m; j < poly->length*m; j++)
         poly->exps[j - i*m] = poly->exps[j];

      _fmpz_mpoly_set_length(poly, poly->length - i, ctx);
   }
}

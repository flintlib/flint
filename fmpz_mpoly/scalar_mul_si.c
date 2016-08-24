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

void _fmpz_mpoly_scalar_mul_si(fmpz * poly1, ulong * exps1,
                    fmpz * poly2, ulong * exps2, slong len2, slong N, slong c)
{
   slong i;

   for (i = 0; i < len2; i++)
      fmpz_mul_si(poly1 + i, poly2 + i, c);

   if (exps1 != exps2)
      mpn_copyi(exps1, exps2, N*len2);
}

void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2,
                                                slong c, fmpz_mpoly_ctx_t ctx)
{
   if (c == 0)
   {
      _fmpz_mpoly_set_length(poly1, 0, ctx);
      return;
   }

   fmpz_mpoly_fit_length(poly1, poly2->length, ctx);

   _fmpz_mpoly_scalar_mul_si(poly1->coeffs, poly1->exps, 
                         poly2->coeffs, poly2->exps, poly2->length, ctx->N, c);
      
   _fmpz_mpoly_set_length(poly1, poly2->length, ctx);
}

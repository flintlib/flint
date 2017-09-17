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

void fmpz_mpoly_get_monomial(ulong * exps, const fmpz_mpoly_t poly, 
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
   slong N = words_per_exp(ctx->n, poly->bits);
   int deg, rev;

   degrev_from_ord(deg, rev, ctx->ord);

   mpoly_get_monomial(exps, poly->exps + N*n, poly->bits, ctx->n, deg, rev);
}

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

void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong j;
    ulong * mon;
    TMP_INIT;

    TMP_START;
    fmpz_mpoly_fit_length(poly, 1, ctx);
    fmpz_set_ui(poly->coeffs + 0, 1);

    mon = (ulong *) TMP_ALLOC((ctx->minfo->nvars)*sizeof(ulong));
    for (j = 0; j < ctx->minfo->nvars; j++)
       mon[j] = (j == i);
    mpoly_set_monomial_ui(poly->exps, mon, poly->bits, ctx->minfo);

    TMP_END;
    _fmpz_mpoly_set_length(poly, 1, ctx);
}

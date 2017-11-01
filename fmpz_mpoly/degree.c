/*
    Copyright (C) 2017 Daniel Schultz

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
#include "assert.h"


void fmpz_mpoly_degrees(slong * degs, const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    int deg, rev;
    degrev_from_ord(deg, rev, ctx->ord);
    mpoly_degrees(degs, poly->exps, poly->length, poly->bits, ctx->n, deg, rev);
}

slong fmpz_mpoly_degree(const fmpz_mpoly_t poly, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong * degs, nvars, ret;
    int deg, rev;
    TMP_INIT;

    TMP_START;
    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;
    degs = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    fmpz_mpoly_degrees(degs, poly, ctx);
    ret = degs[var];

    TMP_END;
    return ret;
}


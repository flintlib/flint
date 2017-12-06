/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_gen(nmod_mpoly_t poly, slong i, const nmod_mpoly_ctx_t ctx)
{
    int deg, rev;
    slong j, nvars;
    ulong * mon;
    TMP_INIT;

    degrev_from_ord(deg, rev, ctx->ord);
    nvars = ctx->n - deg;

    nmod_mpoly_fit_length(poly, 1, ctx);

    poly->coeffs[0] = WORD(1);

    TMP_START;

    mon = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    for (j = 0; j < nvars; j++)
       mon[j] = (j == i);
    mpoly_set_monomial(poly->exps, mon, poly->bits, ctx->n, deg, rev);

    TMP_END;

    _nmod_mpoly_set_length(poly, 1, ctx);
}

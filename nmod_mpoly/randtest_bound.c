/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_randtest_bound(nmod_mpoly_t poly, flint_rand_t state,
                     slong length, slong exp_bound, const nmod_mpoly_ctx_t ctx)
{
    slong i, j, nvars = ctx->minfo->nvars;
    ulong * exp;
    TMP_INIT;

    TMP_START;

    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    nmod_mpoly_zero(poly, ctx);
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bound);

        nmod_mpoly_set_term_ui_ui(poly, n_randtest(state), exp, ctx);
    }
}

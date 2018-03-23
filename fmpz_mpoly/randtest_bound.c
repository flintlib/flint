/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_randtest_bound(fmpz_mpoly_t poly, flint_rand_t state,
    slong length, mp_bitcnt_t coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, nvars = ctx->minfo->nvars;
    fmpz_t c;
    ulong * exp;
    TMP_INIT;

    TMP_START;
    fmpz_init(c);

    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    fmpz_mpoly_zero(poly, ctx);
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bound);

        fmpz_randtest(c, state, coeff_bits);

        fmpz_mpoly_set_term_fmpz_ui(poly, c, exp, ctx);
    }    

    fmpz_clear(c);
}

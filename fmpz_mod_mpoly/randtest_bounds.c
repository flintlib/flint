/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A, flint_rand_t state,
             slong length, ulong * exp_bounds, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j, nvars = ctx->minfo->nvars;
    ulong * exp;
    TMP_INIT;

    TMP_START;
    exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    fmpz_mod_mpoly_zero(A, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(A, 0, MPOLY_MIN_BITS, ctx);
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < nvars; j++)
            exp[j] = n_randint(state, exp_bounds[j]);
        _fmpz_mod_mpoly_push_exp_ui(A, exp, ctx);
        fmpz_randm(A->coeffs + A->length - 1, state,
                                            fmpz_mod_ctx_modulus(ctx->ffinfo));
    }

    fmpz_mod_mpoly_sort_terms(A, ctx);
    fmpz_mod_mpoly_combine_like_terms(A, ctx);

    TMP_END;
}

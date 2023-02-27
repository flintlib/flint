/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A, flint_rand_t state,
                                  slong length, flint_bitcnt_t exp_bits,
                                                const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j, nvars = ctx->minfo->nvars;
    fmpz * exp;
    TMP_INIT;

    TMP_START;

    exp = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exp + j);

    fmpz_mod_mpoly_zero(A, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(A, 0, MPOLY_MIN_BITS, ctx);
    for (i = 0; i < length; i++)
    {
        mpoly_monomial_randbits_fmpz(exp, state, exp_bits, ctx->minfo);
        _fmpz_mod_mpoly_push_exp_ffmpz(A, exp, ctx);
        fmpz_randm(A->coeffs + A->length - 1, state,
                                            fmpz_mod_ctx_modulus(ctx->ffinfo));
    }
    fmpz_mod_mpoly_sort_terms(A, ctx);
    fmpz_mod_mpoly_combine_like_terms(A, ctx);

    for (j = 0; j < nvars; j++)
        fmpz_clear(exp + j);

    TMP_END;
}

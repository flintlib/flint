/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

int fmpz_mod_mpoly_equal_fmpz(
    const fmpz_mod_mpoly_t A,
    const fmpz_t c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > 1)
        return 0;

    if (A->length < 1)
        return fmpz_divisible(c, fmpz_mod_mpoly_ctx_modulus(ctx));

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        return 0;

    return fmpz_mod_equal_fmpz(A->coeffs + 0, c, ctx->ffinfo);
}

int fmpz_mod_mpoly_equal_ui(
    const fmpz_mod_mpoly_t A,
    ulong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int result;
    fmpz_t C;
    fmpz_init_set_ui(C, c);
    result = fmpz_mod_mpoly_equal_fmpz(A, C, ctx);
    fmpz_clear(C);
    return result;
}

int fmpz_mod_mpoly_equal_si(
    const fmpz_mod_mpoly_t A,
    slong c,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > 1)
        return 0;

    if (A->length < 1)
    {
        ulong uc;

        if (c == 0)
            return 0;

        if (!fmpz_abs_fits_ui(fmpz_mod_mpoly_ctx_modulus(ctx)))
            return 0;

        uc = FLINT_ABS(c);
        return 0 == (uc % fmpz_get_ui(fmpz_mod_mpoly_ctx_modulus(ctx)));
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        return 0;

    return fmpz_mod_equal_si(A->coeffs + 0, c, ctx->ffinfo);
}

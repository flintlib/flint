/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_monomial_fmpz(fmpz_mpoly_t poly, 
                              slong n, fmpz ** exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, exp_bits, N;
    fmpz * temp_exp;
    TMP_INIT;

    if (n > poly->length)
        flint_throw(FLINT_ERROR, "Invalid index in fmpz_mpoly_set_monomial");

    TMP_START;
    temp_exp = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init_set(temp_exp + i, exp[i]);
    exp_bits = mpoly_exp_bits_required_fmpz(temp_exp, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

    fmpz_mpoly_fit_length(poly, n + 1, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_set_monomial_fmpz(poly->exps + N*n, temp_exp, poly->bits, ctx->minfo);

    _fmpz_mpoly_set_length(poly, FLINT_MAX(n + 1, poly->length), ctx);

    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_clear(temp_exp + i);

    TMP_END;
}

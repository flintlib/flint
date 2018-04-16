/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_termexp_fmpz(fmpz_mpoly_t poly, 
                      slong n, fmpz * const * exps, const fmpz_mpoly_ctx_t ctx)
{
    slong exp_bits, N;

    if ((ulong) n >= (ulong) poly->length)
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_termexp_fmpz");

    exp_bits = mpoly_exp_bits_required_pfmpz(exps, ctx->minfo);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    fmpz_mpoly_fit_bits(poly, exp_bits, ctx);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_set_monomial_pfmpz(poly->exps + N*n, exps, poly->bits, ctx->minfo);
}

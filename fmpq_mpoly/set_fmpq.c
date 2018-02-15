/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_set_fmpz(fmpq_mpoly_t poly,
                                    const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(poly->zpoly->bits, ctx->zctx->minfo);

    if (fmpz_is_zero(c))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    N = mpoly_words_per_exp(poly->zpoly->bits, ctx->zctx->minfo);

    fmpz_mpoly_fit_length(poly->zpoly, 1, ctx->zctx);
    fmpz_set_si(poly->zpoly->coeffs + 0, WORD(1));
    fmpz_set(fmpq_numref(poly->content), c);
    fmpz_set_si(fmpq_denref(poly->content), WORD(1));

    mpoly_monomial_zero(poly->zpoly->exps + N*0, N);
    _fmpz_mpoly_set_length(poly->zpoly, 1, ctx->zctx);
}

void fmpq_mpoly_set_fmpq(fmpq_mpoly_t poly,
                                    const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(poly->zpoly->bits, ctx->zctx->minfo);

    if (fmpq_is_zero(c))
    {
        fmpq_mpoly_zero(poly, ctx);
        return;
    }

    N = mpoly_words_per_exp(poly->zpoly->bits, ctx->zctx->minfo);

    fmpz_mpoly_fit_length(poly->zpoly, 1, ctx->zctx);
    fmpz_set_si(poly->zpoly->coeffs + 0, WORD(1));
    fmpq_set(poly->content, c);

    mpoly_monomial_zero(poly->zpoly->exps + N*0, N);
    _fmpz_mpoly_set_length(poly->zpoly, 1, ctx->zctx);
}

/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                    const fmpz_t k, const fmpq_mpoly_ctx_t ctx)
{
    if (fmpz_sgn(k) < 0)
    {
        flint_throw(FLINT_ERROR, "Negative power in fmpq_mpoly_pow_fmpz");
    }

    /*
        since there is no fmpq_pow_ui function, it is simplest
        to just send the [ui but not si] case down the fmpz path
    */
    if (fmpz_fits_si(k))
    {
        fmpq_pow_si(A->content, B->content, fmpz_get_si(k));
        fmpz_mpoly_pow_ui(A->zpoly, B->zpoly, fmpz_get_si(k), ctx->zctx);
        return;
    }

    /*
        we are raising a polynomial to an unreasonable exponent.
        It must either be zero or a monomial with coefficient +-1.
    */

    if (fmpq_mpoly_is_zero(B, ctx))
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }

    if (B->zpoly->length != WORD(1))
    {
        flint_throw(FLINT_ERROR, "Multinomial in fmpq_mpoly_pow_fmpz");
    }

    if (!fmpq_is_pm1(B->content))
    {
        flint_throw(FLINT_ERROR, "Non-unit coefficient in fmpq_mpoly_pow_fmpz");
    }

    if (fmpq_is_one(B->content) || fmpz_is_even(k))
    {
        fmpq_one(A->content);
    }
    else
    {
        fmpq_one(A->content);
        fmpq_neg(A->content, A->content);
    }

    fmpz_mpoly_pow_fmpz(A->zpoly, B->zpoly, k, ctx->zctx);
}

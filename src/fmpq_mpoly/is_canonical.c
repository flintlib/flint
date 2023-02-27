/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_is_canonical(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t g;
    int ret;

    if (!fmpq_is_canonical(A->content))
        return 0;

    if (!fmpz_mpoly_is_canonical(A->zpoly, ctx->zctx))
        return 0;

    if (fmpq_is_zero(A->content))
        return A->zpoly->length == 0;

    if (A->zpoly->length == 0)
        return !!fmpq_is_zero(A->content);

    if (fmpz_sgn(A->zpoly->coeffs + 0) <= 0)
        return 0;

    fmpz_init(g);
    _fmpz_vec_content(g, A->zpoly->coeffs, A->zpoly->length);
    ret = fmpz_is_one(g);
    fmpz_clear(g);
    return ret;
}


void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t g;

    if (!fmpq_is_canonical(A->content))
        flint_throw(FLINT_ERROR, "Polynomial content is not canonical");

    fmpz_mpoly_assert_canonical(A->zpoly, ctx->zctx);

    if (fmpq_is_zero(A->content))
    {
        if (A->zpoly->length != 0)
            flint_throw(FLINT_ERROR, "Polynomial content is zero but zpoly is not");
        return;
    }

    if (A->zpoly->length == 0)
    {
        if (!fmpq_is_zero(A->content))
            flint_throw(FLINT_ERROR, "Polynomial zpoly is zero but content is not");
        return;
    }

    if (fmpz_sgn(A->zpoly->coeffs + 0) <= 0)
        flint_throw(FLINT_ERROR, "Polynomial zpoly has negative leading coefficient");

    fmpz_init(g);
    _fmpz_vec_content(g, A->zpoly->coeffs, A->zpoly->length);
    if (!fmpz_is_one(g))
        flint_throw(FLINT_ERROR, "Polynomial zpoly has content");

    fmpz_clear(g);
}

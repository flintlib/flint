/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

    if (mpoly_monomials_overflow_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_inorder_test(poly->exps, poly->length, poly->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < poly->length; i++)
    {
        if (fmpz_equal_si(poly->coeffs + i, 0))
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }
}

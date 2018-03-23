/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    const fmpz_mpoly_struct * zpoly = poly->zpoly;
    const fmpz_mpoly_ctx_struct * zctx = ctx->zctx;

    if (mpoly_monomials_overflow_test(zpoly->exps, zpoly->length, zpoly->bits, zctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_inorder_test(zpoly->exps, zpoly->length, zpoly->bits, zctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < zpoly->length; i++)
    {
        if (fmpz_is_zero(zpoly->coeffs + i))
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }

    if (zpoly->length > 0)
    {
        fmpz_t c;

        if (fmpq_is_zero(poly->content))
        {
            flint_throw(FLINT_ERROR, "Polynomial content is zero");
        }

        if (fmpz_sgn(zpoly->coeffs + 0) <= 0)
        {
            flint_throw(FLINT_ERROR, "Polynomial leading term is not stored positive");
        }

        fmpz_init(c);
        _fmpz_vec_content(c, zpoly->coeffs, zpoly->length);
        if (!fmpz_is_pm1(c))
        {
            flint_throw(FLINT_ERROR, "Polynomial stored coefficients are not relatively prime");
        }

        fmpz_clear(c);

    } else {
        if (!fmpq_is_zero(poly->content))
        {
            flint_throw(FLINT_ERROR, "Polynomial is zero but content is not zero");
        }
    }
}

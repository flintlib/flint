/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (fmpz_is_zero(A->coeffs + i))
            return 0;
    }

    return 1;
}

void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    if (mpoly_monomials_overflow_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

    if (!mpoly_monomials_valid_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

    if (!mpoly_monomials_inorder_test(A->exps, A->length, A->bits, ctx->minfo))
        flint_throw(FLINT_ERROR, "Polynomial exponents out of order");

    for (i = 0; i < A->length; i++)
    {
        if (fmpz_is_zero(A->coeffs + i))
            flint_throw(FLINT_ERROR, "Polynomial has a zero coefficient");
    }

    for (i = A->length; i < A->alloc; i++)
    {
        if (COEFF_IS_MPZ(A->coeffs[i]))
            flint_throw(FLINT_ERROR, "Polynomial has a big coeff past length");
    }
}

int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > WORD(1))
        return 0;

    if (A->length == WORD(0))
        return 1;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    return mpoly_monomial_is_zero(A->exps + N*0, N);
}

int fmpz_mpoly_is_gen(const fmpz_mpoly_t A,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    if (A->length != WORD(1))
        return 0;

    if (!fmpz_is_one(A->coeffs + 0))
        return 0;

    return mpoly_is_gen(A->exps, var, A->bits, ctx->minfo);
}

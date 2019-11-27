/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_cmp(const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    if (A->zpoly->length != 1 || B->zpoly->length != 1
        || !fmpq_is_one(A->content) || !fmpq_is_one(B->content))
    {
        flint_throw(FLINT_ERROR, "Inputs to cmp are not both monomials");
    }

    return mpoly_monomial_cmp_general(A->zpoly->exps, A->zpoly->bits,
                                      B->zpoly->exps, B->zpoly->bits,
                                                             ctx->zctx->minfo);
}

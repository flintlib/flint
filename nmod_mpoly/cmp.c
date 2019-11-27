/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    if (A->length != 1 || B->length != 1
        || A->coeffs[0] != UWORD(1) || B->coeffs[0] != UWORD(1))
    {
        flint_throw(FLINT_ERROR, "Inputs to cmp are not both monomials");
    }

    return mpoly_monomial_cmp_general(A->exps, A->bits, B->exps, B->bits,
                                                                   ctx->minfo);
}

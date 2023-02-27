/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"


int fmpz_mod_mpoly_equal(
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    if (A == B)
        return 1;

    if (A->length != B->length)
        return 0;

    if (!_fmpz_vec_equal(A->coeffs, B->coeffs, A->length))
        return 0;

    return 0 == mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits,
                                                        A->length, ctx->minfo);
}

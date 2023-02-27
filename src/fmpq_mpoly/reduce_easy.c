/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/* don't do too much work if length(A) matches easy_length */
void fmpq_mpoly_reduce_easy(
    fmpq_mpoly_t A,
    slong easy_length,
    const fmpq_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(easy_length > 0);
    FLINT_ASSERT(A->zpoly->length <= easy_length);

    if (A->zpoly->length != easy_length)
    {
        fmpq_mpoly_reduce(A, ctx);
    }
    else if (fmpz_sgn(A->zpoly->coeffs + 0) < 0)
    {
        fmpz_neg(fmpq_numref(A->content), fmpq_numref(A->content));
        _fmpz_vec_neg(A->zpoly->coeffs, A->zpoly->coeffs, A->zpoly->length);
    }

    FLINT_ASSERT(fmpq_mpoly_is_canonical(A, ctx));
}

/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_is_monic(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
{
    int res;
    fmpz_t t;

    if (A->zpoly->length < 1)
        return 0;

    if (fmpz_is_one(fmpq_numref(A->content)) &&
        fmpz_equal(fmpq_denref(A->content), A->zpoly->coeffs + 0))
    {
        return 1;
    }

    fmpz_init(t);
    fmpz_mul(t, fmpq_numref(A->content), A->zpoly->coeffs + 0);
    res = fmpz_equal(t, fmpq_denref(A->content));
    fmpz_clear(t);
    return res;
}


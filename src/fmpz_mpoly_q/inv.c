/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(fmpz_mpoly_q_numref(x), ctx))
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_inv: division by zero\n");
    }

    if (res != x)
        fmpz_mpoly_q_set(res, x, ctx);

    fmpz_mpoly_swap(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res), ctx);

    if (fmpz_sgn(fmpz_mpoly_q_denref(res)->coeffs) < 0)
    {
        fmpz_mpoly_neg(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(res), ctx);
        fmpz_mpoly_neg(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(res), ctx);
    }
}


/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

int
fmpz_mpoly_q_is_canonical(const fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
{
    if (!fmpz_mpoly_is_canonical(fmpz_mpoly_q_numref(res), ctx))
        return 0;

    if (!fmpz_mpoly_is_canonical(fmpz_mpoly_q_denref(res), ctx))
        return 0;

    if (fmpz_mpoly_is_zero(fmpz_mpoly_q_denref(res), ctx))
        return 0;

    if (fmpz_sgn(fmpz_mpoly_q_denref(res)->coeffs) < 0)
        return 0;

    {
        int ans;
        fmpz_mpoly_t g;
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_gcd_assert_successful(g, fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res), ctx);

        ans = fmpz_mpoly_is_one(g, ctx);
        fmpz_mpoly_clear(g, ctx);
        return ans;
    }
}


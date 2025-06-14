/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

int
fmpz_mod_mpoly_q_is_canonical(const fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (!fmpz_mod_mpoly_is_canonical(fmpz_mod_mpoly_q_numref(res), ctx))
        return 0;

    if (!fmpz_mod_mpoly_is_canonical(fmpz_mod_mpoly_q_denref(res), ctx))
        return 0;

    {
        int ans;
        int ans_coeff;
        fmpz_mod_mpoly_t g;
        fmpz_mod_mpoly_init(g, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(g, fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res), ctx);

        ans = fmpz_mod_mpoly_is_one(g, ctx);
        ans_coeff = fmpz_mod_is_one(fmpz_mod_mpoly_q_denref(res)->coeffs + 0, ctx->ffinfo);
        fmpz_mod_mpoly_clear(g, ctx);

        return ans && ans_coeff;
    }
}

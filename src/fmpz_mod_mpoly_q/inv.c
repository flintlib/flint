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

void
fmpz_mod_mpoly_q_inv(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (res != x)
        fmpz_mod_mpoly_q_set(res, x, ctx);

    fmpz_mod_mpoly_swap(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res), ctx);

    if (!fmpz_is_one(fmpz_mod_mpoly_q_denref(res)->coeffs))
    {
        fmpz_t g;
        fmpz_init(g);

        fmpz_mod_inv(g, fmpz_mod_mpoly_q_denref(res)->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);
    
        fmpz_clear(g);
    }
}

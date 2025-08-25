/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include "fmpz_mod_mpoly_q.h"

void
fmpz_mod_mpoly_q_canonicalise(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_one(fmpz_mod_mpoly_q_denref(res), ctx))
    {
        return;
    }
    else if (fmpz_mod_mpoly_is_zero(fmpz_mod_mpoly_q_numref(res), ctx))
    {
        fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);
        return;
    }
    else if (fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_denref(res), ctx))
    {
        fmpz_t g;
        fmpz_init(g);

        fmpz_mod_inv(g, fmpz_mod_mpoly_q_denref(res)->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
        fmpz_mod_mpoly_one(fmpz_mod_mpoly_q_denref(res), ctx);

        fmpz_clear(g);
        return;
    }
    else if (fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_numref(res), ctx))
    {
        fmpz_t g;
        fmpz_init(g);

        fmpz_mod_inv(g, fmpz_mod_mpoly_q_denref(res)->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);

        fmpz_clear(g);
        return;
    }
    else
    {
        fmpz_t g;
        fmpz_mod_mpoly_t t;

        fmpz_init(g);
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(t, fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res), ctx);

        _fmpz_mod_mpoly_q_mpoly_divexact(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), t, ctx);
        _fmpz_mod_mpoly_q_mpoly_divexact(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), t, ctx);

        if (!fmpz_is_one(fmpz_mod_mpoly_q_denref(res)->coeffs))
        {
            fmpz_mod_inv(g, fmpz_mod_mpoly_q_denref(res)->coeffs, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);
        }

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_clear(g);
    }
}


/*
    Copyright (C) 2020 Fredrik Johansson

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

        _fmpz_vec_content(g, fmpz_mod_mpoly_q_numref(res)->coeffs, fmpz_mod_mpoly_q_numref(res)->length);
        fmpz_gcd(g, g, fmpz_mod_mpoly_q_denref(res)->coeffs);

        if (fmpz_sgn(fmpz_mod_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_neg(g, g);

        if (!fmpz_is_one(g))
        {
            fmpz_mod_inv(g, g, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);
        }

        fmpz_clear(g);
    }
    else if (fmpz_mod_mpoly_is_fmpz(fmpz_mod_mpoly_q_numref(res), ctx))
    {
        fmpz_t g;
        fmpz_init(g);

        _fmpz_vec_content(g, fmpz_mod_mpoly_q_denref(res)->coeffs, fmpz_mod_mpoly_q_denref(res)->length);
        fmpz_gcd(g, g, fmpz_mod_mpoly_q_numref(res)->coeffs);

        if (fmpz_sgn(fmpz_mod_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_neg(g, g);

        if (!fmpz_is_one(g))
        {
            fmpz_mod_inv(g, g, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
        
        }

        fmpz_clear(g);
    }
    else
    {
        fmpz_t n, d, g_coeff;
        fmpz_mod_mpoly_t g;
        fmpz_init(d);
        fmpz_init(n);
        fmpz_init(g_coeff);
        fmpz_mod_mpoly_init(g, ctx);

        _fmpz_vec_content(n, fmpz_mod_mpoly_q_numref(res)->coeffs, fmpz_mod_mpoly_q_numref(res)->length);
        _fmpz_vec_content(d, fmpz_mod_mpoly_q_denref(res)->coeffs, fmpz_mod_mpoly_q_denref(res)->length);
        fmpz_gcd(g_coeff, n, d);

        fmpz_mod_mpoly_gcd_assert_successful(g, fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res), ctx);

        if (fmpz_sgn(fmpz_mod_mpoly_q_denref(res)->coeffs) < 0)
            fmpz_mod_mpoly_neg(g, g, ctx);

        if (!fmpz_mod_mpoly_is_one(g, ctx))
        {
            fmpz_mod_inv(g_coeff, g_coeff, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g_coeff, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g_coeff, ctx);
            fmpz_mod_mpoly_divides(fmpz_mod_mpoly_q_denref(res), fmpz_mod_mpoly_q_denref(res), g, ctx);
            fmpz_mod_mpoly_divides(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_numref(res), g, ctx);
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_clear(n);
        fmpz_clear(d);
        fmpz_clear(g_coeff);
    }
}

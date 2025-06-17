/*
    Copyright (C) 2020 Fredrik Johansson
    Copyright (C) 2025 Andrii Yanovets

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mod_mpoly_q.h"

void
_fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(y_num, ctx))
    {
        flint_throw(FLINT_ERROR, "_fmpz_mod_mpoly_q_div: division by zero\n");
    }

    if (fmpz_mod_mpoly_is_zero(x_num, ctx) || fmpz_mod_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mod_mpoly_zero(res_num, ctx);
        fmpz_mod_mpoly_one(res_den, ctx);
        return;
    }

    if (res_num == y_num)
    {
        fmpz_mod_mpoly_t t, u;
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        _fmpz_mod_mpoly_q_mul(t, u, x_num, x_den, y_den, y_num, ctx);
        fmpz_mod_mpoly_swap(res_num, t, ctx);
        fmpz_mod_mpoly_swap(res_den, u, ctx);

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
    }
    else
    {
        _fmpz_mod_mpoly_q_mul(res_num, res_den, x_num, x_den, y_den, y_num, ctx);
    }

    if (!fmpz_is_one(res_den->coeffs))
    {
        fmpz_t g;
        fmpz_init(g);

        fmpz_mod_inv(g, res_den->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, g, ctx);
        fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_den, res_den, g, ctx);
    
        fmpz_clear(g);
    }
}

void
fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpz_mod_mpoly_q_numref(y), fmpz_mod_mpoly_q_denref(y),
                ctx);
}

int
fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    int invertible;
    fmpz_init(t);

    fmpz_mod_set_fmpz(t, fmpq_numref(y), ctx->ffinfo);
    invertible = !fmpz_is_zero(t);

    if (invertible)
    {
        fmpz_mod_inv(t, t, ctx->ffinfo);
        fmpz_mod_mul_fmpz(t, t, fmpq_denref(y), ctx->ffinfo);
        fmpz_mod_mpoly_q_mul_fmpz_mod(res, x, t, ctx);
    }

    fmpz_clear(t);
    return invertible;
}

int
fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    int invertible;
    fmpz_init(t);

    fmpz_mod_set_fmpz(t, y, ctx->ffinfo);
    invertible = !fmpz_is_zero(t);

    if (invertible)
    {
        fmpz_mod_inv(t, t, ctx->ffinfo);
        fmpz_mod_mpoly_q_mul_fmpz_mod(res, x, t, ctx);
    }

    fmpz_clear(t);
    return invertible;
}


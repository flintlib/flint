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
_fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mod_mpoly_neg(res_num, y_num, ctx);
        fmpz_mod_mpoly_set(res_den, y_den, ctx);
        return;
    }

    if (fmpz_mod_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mod_mpoly_set(res_num, x_num, ctx);
        fmpz_mod_mpoly_set(res_den, x_den, ctx);
        return;
    }

    if (fmpz_mod_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mod_mpoly_sub(res_num, x_num, y_num, ctx);

        if (fmpz_mod_mpoly_is_one(x_den, ctx) || fmpz_mod_mpoly_is_zero(res_num, ctx))
        {
            fmpz_mod_mpoly_one(res_den, ctx);
            return;
        }
        else
        {
            fmpz_mod_mpoly_t t;
            fmpz_mod_mpoly_init(t, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, x_den, ctx);

            if (fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                _fmpz_mod_mpoly_q_mpoly_divexact(res_num, res_num, t, ctx);
                _fmpz_mod_mpoly_q_mpoly_divexact(res_den, x_den, t, ctx);
            }

            fmpz_mod_mpoly_clear(t, ctx);
            return;
        }
    }

    if (fmpz_mod_mpoly_is_one(x_den, ctx))
    {  
        if (res_num == y_num)
        {
            fmpz_mod_mpoly_t t;
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mod_mpoly_sub(res_num, t, y_num, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mod_mpoly_mul(res_num, x_num, y_den, ctx);
            fmpz_mod_mpoly_sub(res_num, res_num, y_num, ctx);
        }
        fmpz_mod_mpoly_set(res_den, y_den, ctx);
        return;
    }

    if (fmpz_mod_mpoly_is_one(y_den, ctx))
    {
        if (res_num == x_num)
        {
            fmpz_mod_mpoly_t t;
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_mul(t, y_num, x_den, ctx);
            fmpz_mod_mpoly_sub(res_num, x_num, t, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mod_mpoly_mul(res_num, y_num, x_den, ctx);
            fmpz_mod_mpoly_sub(res_num, x_num, res_num, ctx);
        }
        fmpz_mod_mpoly_set(res_den, x_den, ctx);
        return;
    }

    {
        fmpz_mod_mpoly_t g, t, u;

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(g, x_den, y_den, ctx);

        if (fmpz_mod_mpoly_is_one(g, ctx))
        {   
            fmpz_mod_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mod_mpoly_mul(u, y_num, x_den, ctx);
            fmpz_mod_mpoly_sub(res_num, t, u, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mod_mpoly_t a, b;
            fmpz_mod_mpoly_init(a, ctx);
            fmpz_mod_mpoly_init(b, ctx);

            _fmpz_mod_mpoly_q_mpoly_divexact(a, x_den, g, ctx);
            _fmpz_mod_mpoly_q_mpoly_divexact(b, y_den, g, ctx);

            fmpz_mod_mpoly_mul(t, x_num, b, ctx);
            fmpz_mod_mpoly_mul(u, y_num, a, ctx);
            fmpz_mod_mpoly_sub(res_num, t, u, ctx);

            fmpz_mod_mpoly_gcd_assert_successful(t, res_num, g, ctx);

            if (fmpz_mod_mpoly_is_one(t, ctx))
            {
                fmpz_mod_mpoly_mul(res_den, x_den, b, ctx);
            }
            else
            {
                _fmpz_mod_mpoly_q_mpoly_divexact(res_num, res_num, t, ctx);
                _fmpz_mod_mpoly_q_mpoly_divexact(g, x_den, t, ctx);
                fmpz_mod_mpoly_mul(res_den, g, b, ctx);
            }

            fmpz_mod_mpoly_clear(a, ctx);
            fmpz_mod_mpoly_clear(b, ctx);
        }

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        return;
    }
}

void
fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpz_mod_mpoly_q_numref(y), fmpz_mod_mpoly_q_denref(y),
                ctx);
}

void
_fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_t y,
            const fmpz_mod_mpoly_ctx_t ctx);

void
fmpz_mod_mpoly_q_sub_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_neg(t, y, ctx->ffinfo);
    _fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                t, ctx);
    fmpz_clear(t);
}

void
fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_set_fmpz(t, y, ctx->ffinfo);
    fmpz_mod_neg(t, t, ctx->ffinfo);
    fmpz_mod_mpoly_q_add_fmpz_mod(res, x, t, ctx);
    fmpz_clear(t);
}

int
fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_is_one(fmpq_denref(y)))
    {
        fmpz_mod_mpoly_q_sub_fmpz(res, x, fmpq_numref(y), ctx);
        return 1;
    }
    else
    {
        fmpz_t t;
        int invertible;
        fmpz_init(t);

        fmpz_mod_set_fmpz(t, fmpq_denref(y), ctx->ffinfo);
        invertible = !fmpz_is_zero(t);

        if (invertible)
        {
            fmpz_mod_inv(t, t, ctx->ffinfo);
            fmpz_mod_mul_fmpz(t, t, fmpq_numref(y), ctx->ffinfo);
            fmpz_mod_mpoly_q_sub_fmpz_mod(res, x, t, ctx);
        }

        fmpz_clear(t);
        return invertible;
    }
}


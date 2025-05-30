/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

void
_fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    if (fmpz_mod_mpoly_is_zero(x_num, ctx) || fmpz_mod_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mod_mpoly_zero(res_num, ctx);
        fmpz_mod_mpoly_one(res_den, ctx);
        return;
    }

    if (fmpz_mod_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mod_mpoly_mul(res_num, x_num, y_num, ctx);
        fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
        return;
    }

    /* todo: special-case integer denominators; scalar multiplication */

    if (fmpz_mod_mpoly_is_one(x_den, ctx))
    {
        fmpz_mod_mpoly_t t;
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(t, x_num, y_den, ctx);

        if (fmpz_mod_mpoly_is_one(t, ctx))
        {
            fmpz_mod_mpoly_mul(res_num, x_num, y_num, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mod_mpoly_t u;
            fmpz_mod_mpoly_init(u, ctx);

            _fmpz_mod_mpoly_q_mpoly_divexact(u, x_num, t, ctx);
            fmpz_mod_mpoly_mul(res_num, u, y_num, ctx);
            _fmpz_mod_mpoly_q_mpoly_divexact(u, y_den, t, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, u, ctx);

            fmpz_mod_mpoly_clear(u, ctx);
        }

        fmpz_mod_mpoly_clear(t, ctx);
        return;
    }

    if (fmpz_mod_mpoly_is_one(y_den, ctx))
    {
        fmpz_mod_mpoly_t t;
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(t, y_num, x_den, ctx);

        if (fmpz_mod_mpoly_is_one(t, ctx))
        {
            fmpz_mod_mpoly_mul(res_num, x_num, y_num, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mod_mpoly_t u;
            fmpz_mod_mpoly_init(u, ctx);

            _fmpz_mod_mpoly_q_mpoly_divexact(u, y_num, t, ctx);
            fmpz_mod_mpoly_mul(res_num, u, x_num, ctx);
            _fmpz_mod_mpoly_q_mpoly_divexact(u, x_den, t, ctx);
            fmpz_mod_mpoly_mul(res_den, y_den, u, ctx);
            
            fmpz_mod_mpoly_clear(u, ctx);
        }

        fmpz_mod_mpoly_clear(t, ctx);
        
        return;
    }

    {
        fmpz_mod_mpoly_t t, u, x, y;

        fmpz_mod_mpoly_init(t, ctx);
        fmpz_mod_mpoly_init(u, ctx);
        fmpz_mod_mpoly_init(x, ctx);
        fmpz_mod_mpoly_init(y, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(t, x_num, y_den, ctx);

        if (fmpz_mod_mpoly_is_one(t, ctx))
        {
            fmpz_mod_mpoly_gcd_assert_successful(u, x_den, y_num, ctx);

            if (fmpz_mod_mpoly_is_one(u, ctx))
            {
                fmpz_mod_mpoly_mul(res_num, x_num, y_num, ctx);
                fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);
            }
            else
            {
                fmpz_mod_mpoly_div(y, y_num, u, ctx);
                fmpz_mod_mpoly_mul(res_num, x_num, y, ctx);

                fmpz_mod_mpoly_div(x, x_den, u, ctx);
                fmpz_mod_mpoly_mul(res_den, x, y_den, ctx);
            }
        }
        else
        {
            fmpz_mod_mpoly_gcd_assert_successful(u, x_den, y_num, ctx);

            if (fmpz_mod_mpoly_is_one(u, ctx))
            {
                fmpz_mod_mpoly_div(x, x_num, t, ctx);
                fmpz_mod_mpoly_mul(res_num, x, y_num, ctx);

                fmpz_mod_mpoly_div(y, y_den, t, ctx);
                fmpz_mod_mpoly_mul(res_den, x_den, y, ctx);
            }
            else
            {
                fmpz_mod_mpoly_div(x, x_num, t, ctx);
                fmpz_mod_mpoly_div(y, y_num, u, ctx);
                fmpz_mod_mpoly_mul(res_num, x, y, ctx);

                fmpz_mod_mpoly_div(x, x_den, u, ctx);
                fmpz_mod_mpoly_div(y, y_den, t, ctx);
                fmpz_mod_mpoly_mul(res_den, x, y, ctx);
            }
        }

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(u, ctx);
        fmpz_mod_mpoly_clear(x, ctx);
        fmpz_mod_mpoly_clear(y, ctx);
    }
}

void
_fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{   
    fmpz_t yy, yy_num, yy_den;
    fmpz_init(yy_num);
    fmpz_init(yy_den);
    fmpz_init(yy);
    fmpz_mod_set_fmpz(yy_den, y_den, ctx->ffinfo);
    fmpz_mod_inv(yy_den, yy_den, ctx->ffinfo);
    fmpz_mod_set_fmpz(yy_num, y_num, ctx->ffinfo);
    fmpz_mod_mul(yy, yy_num, yy_den, ctx->ffinfo);

    if (fmpz_mod_mpoly_is_zero(x_num, ctx) || fmpz_is_zero(yy))
    {
        fmpz_mod_mpoly_zero(res_num, ctx);
        fmpz_mod_mpoly_one(res_den, ctx);

        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
        return;
    }
    
    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, x_num, yy, ctx);
    fmpz_mod_mpoly_set(res_den, x_den, ctx);

    fmpz_clear(yy_num);
    fmpz_clear(yy_den);
    fmpz_clear(yy);
    return;    
    
}

void
fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpz_mod_mpoly_q_numref(y), fmpz_mod_mpoly_q_denref(y),
                ctx);
}

void
fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}

void
fmpz_mod_mpoly_q_mul_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t one;
    *one = 1;
    _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                y, one,
                ctx);
}

/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_q.h"

static void
_fmpz_mod_mpoly_q_sub_fmpz_mod_mpoly_den(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
            const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den,
            const fmpz_mod_mpoly_t y_num, const fmpz_t y_den,
            const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t g;
    fmpz_init(g);

    if (fmpz_mod_mpoly_is_fmpz(y_num, ctx))
    {
        if (res_num == x_num || res_num == y_num)
        {
            fmpz_t t, u;
            fmpz_init_set(t, y_num->coeffs);
            fmpz_init_set(u, y_den);
            _fmpz_mod_mpoly_q_sub_fmpq(res_num, res_den, x_num, x_den, t, u, ctx);
            fmpz_clear(t);
            fmpz_clear(u);
        }
        else
        {
            _fmpz_mod_mpoly_q_sub_fmpq(res_num, res_den, x_num, x_den, y_num->coeffs, y_den, ctx);
        }
        return;
    }

    if (fmpz_mod_mpoly_is_fmpz(x_den, ctx))
    {
        fmpz_gcd(g, x_den->coeffs, y_den);

        if (fmpz_mod_is_one(g, ctx->ffinfo))
        {
            fmpz_mod_mpoly_t t, u;

            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_sub does not work inplace */
            fmpz_mod_mpoly_scalar_mul_fmpz(t, y_num, x_den->coeffs, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(u, x_num, y_den, ctx);
            fmpz_mod_mpoly_sub(res_num, u, t, ctx);
            fmpz_mod_mul(g, x_den->coeffs, y_den, ctx->ffinfo);
            fmpz_mod_mpoly_set_fmpz(res_den, g, ctx);

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_t a, b, g_inv;
            fmpz_mod_mpoly_t t, u;

            fmpz_init(a);
            fmpz_init(b);
            fmpz_init(g_inv);
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            fmpz_mod_inv(g_inv, g, ctx->ffinfo);
            fmpz_mod_mul(a, y_den, g_inv, ctx->ffinfo);
            fmpz_mod_mul(b, x_den->coeffs, g_inv, ctx->ffinfo);

            fmpz_mod_mpoly_scalar_mul_fmpz(t, y_num, b, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(u, x_num, a, ctx);
            fmpz_mod_mpoly_sub(res_num, u, t, ctx);

            if (fmpz_mod_mpoly_is_zero(res_num, ctx))
                fmpz_one(a);
            else
                _fmpz_vec_content2(a, res_num->coeffs, res_num->length, g);

            if (fmpz_mod_is_one(a, ctx->ffinfo))
            {
                fmpz_mod_mul(g, b, y_den, ctx->ffinfo);
                fmpz_mod_mpoly_set_fmpz(res_den, g, ctx);
            }
            else
            {   
                fmpz_mod_inv(a, a, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, a, ctx);
                fmpz_mod_mul(g, y_den, a, ctx->ffinfo);
                fmpz_mod_mul(g, g, b, ctx->ffinfo);
                fmpz_mod_mpoly_set_fmpz(res_den, g, ctx);
            }

            fmpz_clear(a);
            fmpz_clear(b);
            fmpz_clear(g_inv);
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
    }
    else
    {
        _fmpz_vec_content2(g, x_den->coeffs, x_den->length, y_den);

        if (fmpz_mod_is_one(g, ctx->ffinfo))
        {
            fmpz_mod_mpoly_t t, u;

            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_sub does not work inplace */
            fmpz_mod_mpoly_mul(t, y_num, x_den, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(u, x_num, y_den, ctx);
            fmpz_mod_mpoly_sub(res_num, u, t, ctx);
            fmpz_mod_set_fmpz(g, y_den, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz(res_den, x_den, g, ctx);

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_t a, g_inv;
            fmpz_mod_mpoly_t b, t, u;

            fmpz_init(a);
            fmpz_init(g_inv);
            fmpz_mod_mpoly_init(b, ctx);
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            fmpz_mod_inv(g, g_inv, ctx->ffinfo);
            fmpz_mod_mul(a, y_den, g, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(b, x_den, g_inv, ctx);

            fmpz_mod_mpoly_mul(t, y_num, b, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(u, x_num, a, ctx);
            fmpz_mod_mpoly_sub(res_num, u, t, ctx);

            if (fmpz_mod_mpoly_is_zero(res_num, ctx))
                fmpz_one(a);
            else
                _fmpz_vec_content2(a, res_num->coeffs, res_num->length, g);

            if (fmpz_mod_is_one(a, ctx->ffinfo))
            {
                fmpz_set(g, y_den);
                fmpz_mod_mpoly_scalar_mul_fmpz(res_den, b, g, ctx);
            }
            else
            {
                fmpz_mod_inv(a, a, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, a, ctx);
                fmpz_mod_mul(g, y_den, a, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz(res_den, b, g, ctx);
            }

            fmpz_clear(a);
            fmpz_clear(g_inv);
            fmpz_mod_mpoly_clear(b, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
    }

    fmpz_clear(g);
}

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
        }
        else if (fmpz_mod_mpoly_is_fmpz(x_den, ctx))
        {
            fmpz_t t;
            fmpz_init(t);

            _fmpz_vec_content2(t, res_num->coeffs, res_num->length, x_den->coeffs);

            if (fmpz_mod_is_one(t, ctx->ffinfo))
            {
                fmpz_mod_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                fmpz_mod_inv(t, t, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, t, ctx);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_den, x_den, t, ctx);
            }

            fmpz_clear(t);
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
        }

        return;
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

    if (fmpz_mod_mpoly_is_fmpz(y_den, ctx))
    {
        _fmpz_mod_mpoly_q_sub_fmpz_mod_mpoly_den(res_num, res_den, x_num, x_den, y_num, y_den->coeffs, ctx);
        return;
    }

    if (fmpz_mod_mpoly_is_fmpz(x_den, ctx))
    {
        _fmpz_mod_mpoly_q_sub_fmpz_mod_mpoly_den(res_num, res_den, y_num, y_den, x_num, x_den->coeffs, ctx);
        fmpz_mod_mpoly_neg(res_num, res_num, ctx);
        return;
    }

    {
        fmpz_mod_mpoly_t g;
        fmpz_mod_mpoly_init(g, ctx);

        fmpz_mod_mpoly_gcd_assert_successful(g, x_den, y_den, ctx);

        if (fmpz_mod_mpoly_is_one(g, ctx))
        {
            fmpz_mod_mpoly_t t, u;

            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_sub does not work inplace */
            fmpz_mod_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mod_mpoly_mul(u, y_num, x_den, ctx);
            fmpz_mod_mpoly_sub(res_num, t, u, ctx);
            fmpz_mod_mpoly_mul(res_den, x_den, y_den, ctx);

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_mod_mpoly_t a, b, t, u;

            fmpz_mod_mpoly_init(a, ctx);
            fmpz_mod_mpoly_init(b, ctx);
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

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
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }

        fmpz_mod_mpoly_clear(g, ctx);
    }
}

void
_fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den,
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

    if (fmpz_is_zero(yy))
    {
        fmpz_mod_mpoly_set(res_num, x_num, ctx);
        fmpz_mod_mpoly_set(res_den, x_den, ctx);
        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
        return;
    }
    fmpz_clear(yy_num);
    fmpz_clear(yy_den);
    fmpz_clear(yy);

    if (fmpz_mod_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mod_mpoly_set_fmpz(res_num, y_num, ctx);
        fmpz_mod_neg(res_num->coeffs, res_num->coeffs, ctx->ffinfo);
        fmpz_mod_mpoly_set_fmpz(res_den, y_den, ctx);
        return;
    }

    /* todo: special-case integer x_den */

    if (fmpz_mod_mpoly_equal_fmpz(x_den, y_den, ctx))
    {
        fmpz_mod_mpoly_sub_fmpz(res_num, x_num, y_num, ctx);

        if (fmpz_mod_is_one(y_den, ctx->ffinfo))
        {
            fmpz_mod_mpoly_one(res_den, ctx);
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);

            _fmpz_vec_content2(t, res_num->coeffs, res_num->length, y_den);

            if (fmpz_mod_is_one(t, ctx->ffinfo))
            {
                fmpz_mod_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                fmpz_mod_inv(t, t, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, t, ctx);
                fmpz_mod_mul(t, y_den, t, ctx->ffinfo);
                fmpz_mod_mpoly_set_fmpz(res_den, t, ctx);
            }

            fmpz_clear(t);
        }

        return;
    }

    if (fmpz_mod_mpoly_is_one(x_den, ctx))
    {
        fmpz_mod_mpoly_scalar_mul_fmpz(res_num, x_num, y_den, ctx);
        fmpz_mod_mpoly_sub_fmpz(res_num, res_num, y_num, ctx);
        fmpz_mod_mpoly_set_fmpz(res_den, y_den, ctx);
        return;
    }

    if (fmpz_mod_is_one(y_den, ctx->ffinfo))
    {
        if (res_num == x_num)
        {
            fmpz_mod_mpoly_t t;
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(t, x_den, y_num, ctx);
            fmpz_mod_mpoly_sub(res_num, x_num, t, ctx);
            fmpz_mod_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mod_mpoly_scalar_mul_fmpz(res_num, x_den, y_num, ctx);
            fmpz_mod_mpoly_sub(res_num, x_num, res_num, ctx);
        }
        fmpz_mod_mpoly_set(res_den, x_den, ctx);
        return;
    }

    {
        fmpz_t g;
        fmpz_init(g);

        _fmpz_vec_content2(g, x_den->coeffs, x_den->length, y_den);

        if (fmpz_mod_is_one(g, ctx->ffinfo))
        {
            fmpz_mod_mpoly_t t, u;

            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            fmpz_mod_mpoly_scalar_mul_fmpz(t, x_num, y_den, ctx);  /* todo: avoid one alloc? */
            fmpz_mod_mpoly_scalar_mul_fmpz(u, x_den, y_num, ctx);
            fmpz_mod_mpoly_sub(res_num, t, u, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);

            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_mod_mpoly_t t, u;
            fmpz_t b, c, g_inv;

            fmpz_init(b);
            fmpz_init(c);
            fmpz_init(g_inv);
            fmpz_mod_mpoly_init(t, ctx);
            fmpz_mod_mpoly_init(u, ctx);

            fmpz_mod_inv(g_inv, g_inv, ctx->ffinfo);
            fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(u, x_den, g_inv, ctx);
            fmpz_mod_mul(b, y_den, g_inv, ctx->ffinfo);

            fmpz_mod_mpoly_scalar_mul_fmpz(t, x_num, b, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(u, u, y_num, ctx);
            fmpz_mod_mpoly_sub(res_num, t, u, ctx);

            _fmpz_vec_content2(c, res_num->coeffs, res_num->length, g);

            if (fmpz_mod_is_one(c, ctx->ffinfo))
            {
                fmpz_mod_mpoly_scalar_mul_fmpz(res_den, x_den, b, ctx);
            }
            else
            {
                fmpz_mod_inv(c, c, ctx->ffinfo);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_num, res_num, c, ctx);
                fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(res_den, x_den, c, ctx);
                fmpz_mod_mpoly_scalar_mul_fmpz(res_den, res_den, b, ctx);
            }

            fmpz_clear(b);
            fmpz_clear(c);
            fmpz_clear(g_inv);
            fmpz_mod_mpoly_clear(t, ctx);
            fmpz_mod_mpoly_clear(u, ctx);
        }

        fmpz_clear(g);
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
fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t one;
    *one = 1;
    _fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                y, one,
                ctx);
}

void
fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    _fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}

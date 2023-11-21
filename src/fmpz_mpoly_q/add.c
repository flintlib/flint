/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

static void
_fmpz_mpoly_q_add_fmpz_mpoly_den(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t g;
    fmpz_init(g);

    if (fmpz_mpoly_is_fmpz(y_num, ctx))
    {
        if (res_num == x_num || res_num == y_num)
        {
            fmpz_t t, u;
            fmpz_init_set(t, y_num->coeffs);
            fmpz_init_set(u, y_den);
            _fmpz_mpoly_q_add_fmpq(res_num, res_den, x_num, x_den, t, u, ctx);
            fmpz_clear(t);
            fmpz_clear(u);
        }
        else
        {
            _fmpz_mpoly_q_add_fmpq(res_num, res_den, x_num, x_den, y_num->coeffs, y_den, ctx);
        }
        return;
    }

    if (fmpz_mpoly_is_fmpz(x_den, ctx))
    {
        fmpz_gcd(g, x_den->coeffs, y_den);

        if (fmpz_is_one(g))
        {
            fmpz_mpoly_t t, u;

            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_add does not work inplace */
            fmpz_mpoly_scalar_mul_fmpz(t, y_num, x_den->coeffs, ctx);
            fmpz_mpoly_scalar_mul_fmpz(u, x_num, y_den, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);
            fmpz_mul(g, x_den->coeffs, y_den);
            fmpz_mpoly_set_fmpz(res_den, g, ctx);

            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_t a, b;
            fmpz_mpoly_t t, u;

            fmpz_init(a);
            fmpz_init(b);
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            fmpz_divexact(a, y_den, g);
            fmpz_divexact(b, x_den->coeffs, g);

            fmpz_mpoly_scalar_mul_fmpz(t, y_num, b, ctx);
            fmpz_mpoly_scalar_mul_fmpz(u, x_num, a, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);

            if (fmpz_mpoly_is_zero(res_num, ctx))
                fmpz_one(a);
            else
                _fmpz_vec_content2(a, res_num->coeffs, res_num->length, g);

            if (fmpz_is_one(a))
            {
                fmpz_mul(g, b, y_den);
                fmpz_mpoly_set_fmpz(res_den, g, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, res_num, a, ctx);
                fmpz_divexact(g, y_den, a);
                fmpz_mul(g, g, b);
                fmpz_mpoly_set_fmpz(res_den, g, ctx);
            }

            fmpz_clear(a);
            fmpz_clear(b);
            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
    }
    else
    {
        _fmpz_vec_content2(g, x_den->coeffs, x_den->length, y_den);

        if (fmpz_is_one(g))
        {
            fmpz_mpoly_t t, u;

            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_add does not work inplace */
            fmpz_mpoly_mul(t, y_num, x_den, ctx);
            fmpz_mpoly_scalar_mul_fmpz(u, x_num, y_den, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);
            fmpz_set(g, y_den);
            fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, g, ctx);

            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_t a;
            fmpz_mpoly_t b, t, u;

            fmpz_init(a);
            fmpz_mpoly_init(b, ctx);
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            fmpz_divexact(a, y_den, g);
            fmpz_mpoly_scalar_divexact_fmpz(b, x_den, g, ctx);

            fmpz_mpoly_mul(t, y_num, b, ctx);
            fmpz_mpoly_scalar_mul_fmpz(u, x_num, a, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);

            if (fmpz_mpoly_is_zero(res_num, ctx))
                fmpz_one(a);
            else
                _fmpz_vec_content2(a, res_num->coeffs, res_num->length, g);

            if (fmpz_is_one(a))
            {
                fmpz_set(g, y_den);
                fmpz_mpoly_scalar_mul_fmpz(res_den, b, g, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, res_num, a, ctx);
                fmpz_divexact(g, y_den, a);
                fmpz_mpoly_scalar_mul_fmpz(res_den, b, g, ctx);
            }

            fmpz_clear(a);
            fmpz_mpoly_clear(b, ctx);
            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
    }

    fmpz_clear(g);
}

void
_fmpz_mpoly_q_add(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mpoly_set(res_num, y_num, ctx);
        fmpz_mpoly_set(res_den, y_den, ctx);
        return;
    }

    if (fmpz_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mpoly_set(res_num, x_num, ctx);
        fmpz_mpoly_set(res_den, x_den, ctx);
        return;
    }

    if (fmpz_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mpoly_add(res_num, x_num, y_num, ctx);

        if (fmpz_mpoly_is_one(x_den, ctx) || fmpz_mpoly_is_zero(res_num, ctx))
        {
            fmpz_mpoly_one(res_den, ctx);
        }
        else if (fmpz_mpoly_is_fmpz(x_den, ctx))
        {
            fmpz_t t;
            fmpz_init(t);

            _fmpz_vec_content2(t, res_num->coeffs, res_num->length, x_den->coeffs);

            if (fmpz_is_one(t))
            {
                fmpz_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, res_num, t, ctx);
                fmpz_mpoly_scalar_divexact_fmpz(res_den, x_den, t, ctx);
            }

            fmpz_clear(t);
        }
        else
        {
            fmpz_mpoly_t t;
            fmpz_mpoly_init(t, ctx);

            fmpz_mpoly_gcd_assert_successful(t, res_num, x_den, ctx);

            if (fmpz_mpoly_is_one(t, ctx))
            {
                fmpz_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                _fmpz_mpoly_q_mpoly_divexact(res_num, res_num, t, ctx);
                _fmpz_mpoly_q_mpoly_divexact(res_den, x_den, t, ctx);
            }

            fmpz_mpoly_clear(t, ctx);
        }

        return;
    }

    if (fmpz_mpoly_is_one(x_den, ctx))
    {
        if (res_num == y_num)
        {
            fmpz_mpoly_t t;
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mpoly_add(res_num, t, y_num, ctx);
            fmpz_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mpoly_mul(res_num, x_num, y_den, ctx);
            fmpz_mpoly_add(res_num, res_num, y_num, ctx);
        }
        fmpz_mpoly_set(res_den, y_den, ctx);
        return;
    }

    if (fmpz_mpoly_is_one(y_den, ctx))
    {
        if (res_num == x_num)
        {
            fmpz_mpoly_t t;
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_mul(t, y_num, x_den, ctx);
            fmpz_mpoly_add(res_num, x_num, t, ctx);
            fmpz_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mpoly_mul(res_num, y_num, x_den, ctx);
            fmpz_mpoly_add(res_num, x_num, res_num, ctx);
        }
        fmpz_mpoly_set(res_den, x_den, ctx);
        return;
    }

    if (fmpz_mpoly_is_fmpz(y_den, ctx))
    {
        _fmpz_mpoly_q_add_fmpz_mpoly_den(res_num, res_den, x_num, x_den, y_num, y_den->coeffs, ctx);
        return;
    }

    if (fmpz_mpoly_is_fmpz(x_den, ctx))
    {
        _fmpz_mpoly_q_add_fmpz_mpoly_den(res_num, res_den, y_num, y_den, x_num, x_den->coeffs, ctx);
        return;
    }

    {
        fmpz_mpoly_t g;
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_gcd_assert_successful(g, x_den, y_den, ctx);

        if (fmpz_mpoly_is_one(g, ctx))
        {
            fmpz_mpoly_t t, u;

            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            /* todo: avoid one alloc? not helpful right now because
               fmpz_mpoly_add does not work inplace */
            fmpz_mpoly_mul(t, x_num, y_den, ctx);
            fmpz_mpoly_mul(u, y_num, x_den, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);
            fmpz_mpoly_mul(res_den, x_den, y_den, ctx);

            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_mpoly_t a, b, t, u;

            fmpz_mpoly_init(a, ctx);
            fmpz_mpoly_init(b, ctx);
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            _fmpz_mpoly_q_mpoly_divexact(a, x_den, g, ctx);
            _fmpz_mpoly_q_mpoly_divexact(b, y_den, g, ctx);

            fmpz_mpoly_mul(t, x_num, b, ctx);
            fmpz_mpoly_mul(u, y_num, a, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);

            fmpz_mpoly_gcd_assert_successful(t, res_num, g, ctx);

            if (fmpz_mpoly_is_one(t, ctx))
            {
                fmpz_mpoly_mul(res_den, x_den, b, ctx);
            }
            else
            {
                _fmpz_mpoly_q_mpoly_divexact(res_num, res_num, t, ctx);
                _fmpz_mpoly_q_mpoly_divexact(g, x_den, t, ctx);
                fmpz_mpoly_mul(res_den, g, b, ctx);
            }

            fmpz_mpoly_clear(a, ctx);
            fmpz_mpoly_clear(b, ctx);
            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }

        fmpz_mpoly_clear(g, ctx);
    }
}

void
_fmpz_mpoly_q_add_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mpoly_set_fmpz(res_num, y_num, ctx);
        fmpz_mpoly_set_fmpz(res_den, y_den, ctx);
        return;
    }

    if (fmpz_is_zero(y_num))
    {
        fmpz_mpoly_set(res_num, x_num, ctx);
        fmpz_mpoly_set(res_den, x_den, ctx);
        return;
    }

    /* todo: special-case integer x_den */

    if (fmpz_mpoly_equal_fmpz(x_den, y_den, ctx))
    {
        fmpz_mpoly_add_fmpz(res_num, x_num, y_num, ctx);

        if (fmpz_is_one(y_den))
        {
            fmpz_mpoly_one(res_den, ctx);
        }
        else
        {
            fmpz_t t;
            fmpz_init(t);

            _fmpz_vec_content2(t, res_num->coeffs, res_num->length, y_den);

            if (fmpz_is_one(t))
            {
                fmpz_mpoly_set(res_den, x_den, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, res_num, t, ctx);
                fmpz_divexact(t, y_den, t);
                fmpz_mpoly_set_fmpz(res_den, t, ctx);
            }

            fmpz_clear(t);
        }

        return;
    }

    if (fmpz_mpoly_is_one(x_den, ctx))
    {
        fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, y_den, ctx);
        fmpz_mpoly_add_fmpz(res_num, res_num, y_num, ctx);
        fmpz_mpoly_set_fmpz(res_den, y_den, ctx);
        return;
    }

    if (fmpz_is_one(y_den))
    {
        if (res_num == x_num)
        {
            fmpz_mpoly_t t;
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t, x_den, y_num, ctx);
            fmpz_mpoly_add(res_num, x_num, t, ctx);
            fmpz_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mpoly_scalar_mul_fmpz(res_num, x_den, y_num, ctx);
            fmpz_mpoly_add(res_num, x_num, res_num, ctx);
        }
        fmpz_mpoly_set(res_den, x_den, ctx);
        return;
    }

    {
        fmpz_t g;
        fmpz_init(g);

        _fmpz_vec_content2(g, x_den->coeffs, x_den->length, y_den);

        if (fmpz_is_one(g))
        {
            fmpz_mpoly_t t, u;

            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            fmpz_mpoly_scalar_mul_fmpz(t, x_num, y_den, ctx);  /* todo: avoid one alloc? */
            fmpz_mpoly_scalar_mul_fmpz(u, x_den, y_num, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);
            fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);

            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }
        else
        {
            fmpz_mpoly_t t, u;
            fmpz_t b, c;

            fmpz_init(b);
            fmpz_init(c);
            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            fmpz_mpoly_scalar_divexact_fmpz(u, x_den, g, ctx);
            fmpz_divexact(b, y_den, g);

            fmpz_mpoly_scalar_mul_fmpz(t, x_num, b, ctx);
            fmpz_mpoly_scalar_mul_fmpz(u, u, y_num, ctx);
            fmpz_mpoly_add(res_num, t, u, ctx);

            _fmpz_vec_content2(c, res_num->coeffs, res_num->length, g);

            if (fmpz_is_one(c))
            {
                fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, b, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, res_num, c, ctx);
                fmpz_mpoly_scalar_divexact_fmpz(res_den, x_den, c, ctx);
                fmpz_mpoly_scalar_mul_fmpz(res_den, res_den, b, ctx);
            }

            fmpz_clear(b);
            fmpz_clear(c);
            fmpz_mpoly_clear(t, ctx);
            fmpz_mpoly_clear(u, ctx);
        }

        fmpz_clear(g);
    }
}

void
fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_add(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpz_mpoly_q_numref(y), fmpz_mpoly_q_denref(y),
                ctx);
}

void
fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}

void
fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_t one;
    *one = 1;
    _fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                y, one,
                ctx);
}

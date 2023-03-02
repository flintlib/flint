/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
_fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_is_zero(y_num))
    {
        fmpz_mpoly_set(res_num, x_num, ctx);
        fmpz_mpoly_set(res_den, x_den, ctx);
        return;
    }

    if (fmpz_mpoly_is_zero(x_num, ctx))
    {
        fmpz_mpoly_set_fmpz(res_num, y_num, ctx);
        fmpz_neg(res_num->coeffs, res_num->coeffs);
        fmpz_mpoly_set_fmpz(res_den, y_den, ctx);
        return;
    }

    /* todo: special-case integer x_den */

    if (fmpz_mpoly_equal_fmpz(x_den, y_den, ctx))
    {
        fmpz_mpoly_sub_fmpz(res_num, x_num, y_num, ctx);

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
        fmpz_mpoly_sub_fmpz(res_num, res_num, y_num, ctx);
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
            fmpz_mpoly_sub(res_num, x_num, t, ctx);
            fmpz_mpoly_clear(t, ctx);
        }
        else
        {
            fmpz_mpoly_scalar_mul_fmpz(res_num, x_den, y_num, ctx);
            fmpz_mpoly_sub(res_num, x_num, res_num, ctx);
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
            fmpz_mpoly_sub(res_num, t, u, ctx);
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
            fmpz_mpoly_sub(res_num, t, u, ctx);

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
fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}


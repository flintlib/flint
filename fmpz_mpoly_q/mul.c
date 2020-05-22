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
_fmpz_mpoly_q_mul(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(x_num, ctx) || fmpz_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mpoly_zero(res_num, ctx);
        fmpz_mpoly_one(res_den, ctx);
        return;
    }

    if (fmpz_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mpoly_mul(res_num, x_num, y_num, ctx);
        fmpz_mpoly_mul(res_den, x_den, y_den, ctx);
        return;
    }

    /* todo: special-case integer denominators; scalar multiplication */

    if (fmpz_mpoly_is_one(x_den, ctx))
    {
        fmpz_mpoly_t t, u;
        fmpz_mpoly_init(t, ctx);

        fmpz_mpoly_gcd_assert_successful(t, x_num, y_den, ctx);

        if (fmpz_mpoly_is_one(t, ctx))
        {
            fmpz_mpoly_mul(res_num, x_num, y_num, ctx);
            fmpz_mpoly_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mpoly_init(u, ctx);
            fmpz_mpoly_divexact(u, x_num, t, ctx);
            fmpz_mpoly_mul(res_num, u, y_num, ctx);
            fmpz_mpoly_divexact(u, y_den, t, ctx);
            fmpz_mpoly_mul(res_den, x_den, u, ctx);
            fmpz_mpoly_clear(u, ctx);
        }

        fmpz_mpoly_clear(t, ctx);
        return;
    }

    if (fmpz_mpoly_is_one(y_den, ctx))
    {
        fmpz_mpoly_t t, u;
        fmpz_mpoly_init(t, ctx);

        fmpz_mpoly_gcd_assert_successful(t, y_num, x_den, ctx);

        if (fmpz_mpoly_is_one(t, ctx))
        {
            fmpz_mpoly_mul(res_num, x_num, y_num, ctx);
            fmpz_mpoly_mul(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mpoly_init(u, ctx);
            fmpz_mpoly_divexact(u, y_num, t, ctx);
            fmpz_mpoly_mul(res_num, u, x_num, ctx);
            fmpz_mpoly_divexact(u, x_den, t, ctx);
            fmpz_mpoly_mul(res_den, y_den, u, ctx);
            fmpz_mpoly_clear(u, ctx);
        }

        fmpz_mpoly_clear(t, ctx);
        return;
    }

    {
        fmpz_mpoly_t t, u, x, y;

        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_init(u, ctx);
        fmpz_mpoly_init(x, ctx);
        fmpz_mpoly_init(y, ctx);

        fmpz_mpoly_gcd_assert_successful(t, x_num, y_den, ctx);

        if (fmpz_mpoly_is_one(t, ctx))
        {
            fmpz_mpoly_gcd_assert_successful(u, x_den, y_num, ctx);

            if (fmpz_mpoly_is_one(u, ctx))
            {
                fmpz_mpoly_mul(res_num, x_num, y_num, ctx);
                fmpz_mpoly_mul(res_den, x_den, y_den, ctx);
            }
            else
            {
                fmpz_mpoly_div(y, y_num, u, ctx);
                fmpz_mpoly_mul(res_num, x_num, y, ctx);

                fmpz_mpoly_div(x, x_den, u, ctx);
                fmpz_mpoly_mul(res_den, x, y_den, ctx);
            }
        }
        else
        {
            fmpz_mpoly_gcd_assert_successful(u, x_den, y_num, ctx);

            if (fmpz_mpoly_is_one(u, ctx))
            {
                fmpz_mpoly_div(x, x_num, t, ctx);
                fmpz_mpoly_mul(res_num, x, y_num, ctx);

                fmpz_mpoly_div(y, y_den, t, ctx);
                fmpz_mpoly_mul(res_den, x_den, y, ctx);
            }
            else
            {
                fmpz_mpoly_div(x, x_num, t, ctx);
                fmpz_mpoly_div(y, y_num, u, ctx);
                fmpz_mpoly_mul(res_num, x, y, ctx);

                fmpz_mpoly_div(x, x_den, u, ctx);
                fmpz_mpoly_div(y, y_den, t, ctx);
                fmpz_mpoly_mul(res_den, x, y, ctx);
            }
        }

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(u, ctx);
        fmpz_mpoly_clear(x, ctx);
        fmpz_mpoly_clear(y, ctx);
    }
}

void
fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_mul(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpz_mpoly_q_numref(y), fmpz_mpoly_q_denref(y),
                ctx);
}


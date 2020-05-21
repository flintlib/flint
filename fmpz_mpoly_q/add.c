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

    /* todo: special-case integer denominators; scalar addition */

    if (fmpz_mpoly_equal(x_den, y_den, ctx))
    {
        fmpz_mpoly_add(res_num, x_num, y_num, ctx);

        if (fmpz_mpoly_is_one(x_den, ctx))
        {
            fmpz_mpoly_one(res_den, ctx);
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
                fmpz_mpoly_divexact(res_num, res_num, t, ctx);
                fmpz_mpoly_divexact(res_den, x_den, t, ctx);
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

    {
        fmpz_mpoly_t g;
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_gcd_assert_successful(g, x_den, y_den, ctx);

        if (fmpz_mpoly_is_one(g, ctx))
        {
            fmpz_mpoly_t t, u;

            fmpz_mpoly_init(t, ctx);
            fmpz_mpoly_init(u, ctx);

            fmpz_mpoly_mul(t, x_num, y_den, ctx);  /* todo: avoid one alloc? */
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

            fmpz_mpoly_divexact(a, x_den, g, ctx);
            fmpz_mpoly_divexact(b, y_den, g, ctx);

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
                fmpz_mpoly_divexact(res_num, res_num, t, ctx);
                fmpz_mpoly_divexact(g, x_den, t, ctx);
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
fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_add(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpz_mpoly_q_numref(y), fmpz_mpoly_q_denref(y),
                ctx);
}


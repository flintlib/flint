/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"

void
_fmpz_mpoly_q_div(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(y_num, ctx))
    {
        flint_throw(FLINT_ERROR, "_fmpz_mpoly_q_div: division by zero\n");
    }

    if (fmpz_mpoly_is_zero(x_num, ctx) || fmpz_mpoly_is_zero(y_num, ctx))
    {
        fmpz_mpoly_zero(res_num, ctx);
        fmpz_mpoly_one(res_den, ctx);
        return;
    }

    if (res_num == y_num)
    {
        fmpz_mpoly_t t, u;
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_init(u, ctx);
        _fmpz_mpoly_q_mul(t, u, x_num, x_den, y_den, y_num, ctx);
        fmpz_mpoly_swap(res_num, t, ctx);
        fmpz_mpoly_swap(res_den, u, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(u, ctx);
    }
    else
    {
        _fmpz_mpoly_q_mul(res_num, res_den, x_num, x_den, y_den, y_num, ctx);
    }

    if (fmpz_sgn(res_den->coeffs) < 0)
    {
        fmpz_mpoly_neg(res_num, res_num, ctx);
        fmpz_mpoly_neg(res_den, res_den, ctx);
    }
}

void
fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_div(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpz_mpoly_q_numref(y), fmpz_mpoly_q_denref(y),
                ctx);
}

void
fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpq_is_zero(y))
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_div_fmpq: division by zero\n");
    }
    else
    {
        if (fmpz_sgn(fmpq_numref(y)) > 0)
        {
            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        fmpq_denref(y), fmpq_numref(y),
                        ctx);
        }
        else
        {
            fmpz_t t, u;
            fmpz_init(t);
            fmpz_init(u);
            fmpz_neg(t, fmpq_numref(y));
            fmpz_neg(u, fmpq_denref(y));

            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        u, t,
                        ctx);

            fmpz_clear(t);
            fmpz_clear(u);
        }
    }
}

void
fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_div_fmpz: division by zero\n");
    }
    else
    {
        if (fmpz_sgn(y) > 0)
        {
            fmpz_t one;
            *one = 1;
            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        one, y,
                        ctx);
        }
        else
        {
            fmpz_t t;
            fmpz_t one;
            *one = -1;
            fmpz_init(t);
            fmpz_neg(t, y);

            _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                        fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                        one, t,
                        ctx);

            fmpz_clear(t);
        }
    }
}


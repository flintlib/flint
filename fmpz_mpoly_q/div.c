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
_fmpz_mpoly_q_div(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(y_num, ctx))
    {
        flint_printf("_fmpz_mpoly_q_div: division by zero\n");
        flint_abort();
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


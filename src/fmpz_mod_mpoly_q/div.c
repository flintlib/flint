/*
    Copyright (C) 2020 Fredrik Johansson

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

    if (fmpz_sgn(res_den->coeffs) < 0)
    {
        fmpz_mod_mpoly_neg(res_num, res_num, ctx);
        fmpz_mod_mpoly_neg(res_den, res_den, ctx);
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

void
fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t yy, yy_num, yy_den;
    fmpz_init(yy_num);
    fmpz_init(yy_den);
    fmpz_init(yy);
    fmpz_mod_inv(yy_den, fmpq_denref(y), ctx->ffinfo);
    fmpz_mod_set_fmpz(yy_num, fmpq_numref(y), ctx->ffinfo);
    fmpz_mod_mul(yy, yy_num, yy_den, ctx->ffinfo);

    if (fmpz_is_zero(yy))
    {
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_q_div_fmpq: division by zero\n");
        fmpz_clear(yy_num);
        fmpz_clear(yy_den);
        fmpz_clear(yy);
    }
    else
    {
        if (fmpz_sgn(fmpq_numref(y)) > 0)
        {
            _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                        fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                        fmpq_denref(y), fmpq_numref(y),
                        ctx);
        }
        else
        {
            fmpz_t t, u;
            fmpz_init(t);
            fmpz_init(u);
            fmpz_mod_neg(t, fmpq_numref(y), ctx->ffinfo);
            fmpz_mod_neg(u, fmpq_denref(y), ctx->ffinfo);

            _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                        fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                        u, t,
                        ctx);

            fmpz_clear(t);
            fmpz_clear(u);
        }
    }
    fmpz_clear(yy_num);
    fmpz_clear(yy_den);
    fmpz_clear(yy);
}

void
fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
{
    fmpz_t yy;
    fmpz_init(yy);
    fmpz_mod_set_fmpz(yy, y, ctx->ffinfo);

    if (fmpz_is_zero(yy))
    {
        flint_throw(FLINT_ERROR, "fmpz_mod_mpoly_q_div_fmpz: division by zero\n");
        fmpz_clear(yy);
    }
    else
    {
        if (fmpz_sgn(y) > 0)
        {
            fmpz_t one;
            *one = 1;
            _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                        fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                        one, y,
                        ctx);
        }
        else
        {
            fmpz_t t;
            fmpz_t one;
            *one = -1;
            fmpz_init(t);
            fmpz_mod_neg(t, y, ctx->ffinfo);

            _fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_numref(res), fmpz_mod_mpoly_q_denref(res),
                        fmpz_mod_mpoly_q_numref(x), fmpz_mod_mpoly_q_denref(x),
                        one, t,
                        ctx);

            fmpz_clear(t);
        }
    }
    fmpz_clear(yy);
}

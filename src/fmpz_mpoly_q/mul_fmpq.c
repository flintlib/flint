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
_fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den,
            const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den,
            const fmpz_t y_num, const fmpz_t y_den,
            const fmpz_mpoly_ctx_t ctx)
{
    if (fmpz_mpoly_is_zero(x_num, ctx) || fmpz_is_zero(y_num))
    {
        fmpz_mpoly_zero(res_num, ctx);
        fmpz_mpoly_one(res_den, ctx);
        return;
    }

    if (fmpz_mpoly_equal_fmpz(x_den, y_den, ctx))
    {
        fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, y_num, ctx);
        fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);
        return;
    }

    if (fmpz_mpoly_is_one(x_den, ctx))
    {
        fmpz_t t;
        fmpz_init(t);

        _fmpz_vec_content2(t, x_num->coeffs, x_num->length, y_den);

        if (fmpz_is_one(t))
        {
            fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, y_num, ctx);
            fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mpoly_scalar_divexact_fmpz(res_num, x_num, t, ctx);
            fmpz_mpoly_scalar_mul_fmpz(res_num, res_num, y_num, ctx);
            fmpz_divexact(t, y_den, t);
            fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, t, ctx);
        }

        fmpz_clear(t);
        return;
    }

    if (fmpz_is_one(y_den))
    {
        fmpz_t t;
        fmpz_init(t);

        _fmpz_vec_content2(t, x_den->coeffs, x_den->length, y_num);

        if (fmpz_is_one(t))
        {
            fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, y_num, ctx);
            fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);
        }
        else
        {
            fmpz_mpoly_scalar_divexact_fmpz(res_den, x_den, t, ctx);
            fmpz_mpoly_scalar_mul_fmpz(res_den, res_den, y_den, ctx);
            fmpz_divexact(t, y_num, t);
            fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, t, ctx);
        }

        fmpz_clear(t);
        return;
    }

    {
        fmpz_t t, u;

        fmpz_init(t);
        fmpz_init(u);

        _fmpz_vec_content2(t, x_num->coeffs, x_num->length, y_den);
        _fmpz_vec_content2(u, x_den->coeffs, x_den->length, y_num);

        if (fmpz_is_one(t))
        {
            if (fmpz_is_one(u))
            {
                fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, y_num, ctx);
                fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, y_den, ctx);
            }
            else
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_den, x_den, u, ctx);
                fmpz_mpoly_scalar_mul_fmpz(res_den, res_den, y_den, ctx);
                fmpz_divexact(u, y_num, u);
                fmpz_mpoly_scalar_mul_fmpz(res_num, x_num, u, ctx);
            }
        }
        else
        {
            if (fmpz_is_one(u))
            {
                fmpz_mpoly_scalar_divexact_fmpz(res_num, x_num, t, ctx);
                fmpz_mpoly_scalar_mul_fmpz(res_num, res_num, y_num, ctx);
                fmpz_divexact(t, y_den, t);
                fmpz_mpoly_scalar_mul_fmpz(res_den, x_den, t, ctx);
            }
            else
            {
                fmpz_t v;
                fmpz_init(v);

                fmpz_mpoly_scalar_divexact_fmpz(res_num, x_num, t, ctx);
                fmpz_divexact(v, y_num, u);
                fmpz_mpoly_scalar_mul_fmpz(res_num, res_num, v, ctx);

                fmpz_mpoly_scalar_divexact_fmpz(res_den, x_den, u, ctx);
                fmpz_divexact(v, y_den, t);
                fmpz_mpoly_scalar_mul_fmpz(res_den, res_den, v, ctx);

                fmpz_clear(v);
            }
        }

        fmpz_clear(t);
        fmpz_clear(u);
    }
}

void
fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(res),
                fmpz_mpoly_q_numref(x), fmpz_mpoly_q_denref(x),
                fmpq_numref(y), fmpq_denref(y),
                ctx);
}


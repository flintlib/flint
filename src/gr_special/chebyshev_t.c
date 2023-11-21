/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_special.h"

int
gr_generic_chebyshev_t2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n1bits, i;
    fmpz_t n1;

    /* not implemented */
    if (fmpz_sgn(n) < 0)
        return GR_UNABLE;

    status |= gr_set(a, x, ctx);
    status |= gr_one(b, ctx);

    if (fmpz_sgn(n) < 0)
        return GR_UNABLE;

    if (fmpz_is_zero(n))
    {
        gr_swap(a, b, ctx);
        return status;
    }

    if (fmpz_is_one(n))
        return status;

    fmpz_init(n1);
    fmpz_sub_ui(n1, n, 1);

    n1bits = fmpz_bits(n1);

    for (i = n1bits - 1; i >= 0; i--)
    {
        if (fmpz_tstbit(n1, i))
        {
            status |= gr_mul(b, b, a, ctx);
            status |= gr_mul_two(b, b, ctx);
            status |= gr_sub(b, b, x, ctx);
            status |= gr_mul(a, a, a, ctx);
            status |= gr_mul_two(a, a, ctx);
            status |= gr_sub_ui(a, a, 1, ctx);
        }
        else
        {
            status |= gr_mul(a, a, b, ctx);
            status |= gr_mul_two(a, a, ctx);
            status |= gr_sub(a, a, x, ctx);
            status |= gr_mul(b, b, b, ctx);
            status |= gr_mul_two(b, b, ctx);
            status |= gr_sub_ui(b, b, 1, ctx);
        }
    }

    fmpz_clear(n1);
    return status;
}

int
gr_generic_chebyshev_t_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong nbits, i, r;
    fmpz_t n1;

    if (fmpz_is_zero(n))
        return gr_one(y, ctx);

    if (fmpz_is_one(n))
        return gr_set(y, x, ctx);

    if (fmpz_sgn(n) < 0)
    {
        fmpz_init(n1);
        fmpz_neg(n1, n);
        status = gr_generic_chebyshev_t_fmpz(y, n1, x, ctx);
        fmpz_clear(n1);
        return status;
    }

    if (gr_is_zero(x, ctx) == T_TRUE)
    {
        int c = fmpz_fdiv_ui(n, 4);
        return gr_set_si(y, (c == 0) - (c == 2), ctx);
    }

    if (gr_is_one(x, ctx) == T_TRUE)
        return gr_one(y, ctx);

    if (gr_is_neg_one(x, ctx) == T_TRUE)
        return fmpz_is_even(n) ? gr_one(y, ctx) : gr_neg_one(y, ctx);

    nbits = fmpz_bits(n);
    r = fmpz_val2(n);

    if (nbits == r + 1)
    {
        status |= gr_sqr(y, x, ctx);
        status |= gr_mul_two(y, y, ctx);
        status |= gr_sub_ui(y, y, 1, ctx);
        r -= 1;
    }
    else
    {
        /* we only need one value, so break out final iteration */
        gr_ptr t, u;
        fmpz_t n1;

        GR_TMP_INIT2(t, u, ctx);
        fmpz_init(n1);

        fmpz_tdiv_q_2exp(n1, n, r + 1);
        fmpz_add_ui(n1, n1, 1);

        status |= gr_generic_chebyshev_t2_fmpz(t, u, n1, x, ctx);
        status |= gr_mul(t, t, u, ctx);
        status |= gr_mul_two(t, t, ctx);
        status |= gr_sub(y, t, x, ctx);

        GR_TMP_CLEAR2(t, u, ctx);
        fmpz_clear(n1);
    }

    for (i = 0; i < r; i++)
    {
        status |= gr_sqr(y, y, ctx);
        status |= gr_mul_two(y, y, ctx);
        status |= gr_sub_ui(y, y, 1, ctx);
    }

    return status;
}

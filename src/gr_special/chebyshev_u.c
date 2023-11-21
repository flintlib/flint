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
gr_generic_chebyshev_u2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong nbits, i;
    gr_ptr t, u;

    /* not implemented */
    if (fmpz_sgn(n) < 0)
        return GR_UNABLE;

    if (fmpz_is_zero(n))
    {
        status |= gr_one(a, ctx);
        status |= gr_zero(b, ctx);
        return status;
    }

    status |= gr_mul_two(a, x, ctx);
    status |= gr_one(b, ctx);

    if (fmpz_is_one(n))
        return status;

    nbits = fmpz_bits(n);

    GR_TMP_INIT2(t, u, ctx);

    for (i = nbits - 2; i >= 0; i--)
    {
        status |= gr_add(t, a, b, ctx);
        status |= gr_sub(u, a, b, ctx);

        if (fmpz_tstbit(n, i))
        {
            status |= gr_submul(b, x, a, ctx);
            status |= gr_mul(a, a, b, ctx);
            status |= gr_neg(a, a, ctx);
            status |= gr_mul_two(a, a, ctx);
            status |= gr_mul(b, t, u, ctx);
        }
        else
        {
            status |= gr_submul(a, x, b, ctx);
            status |= gr_mul(b, a, b, ctx);
            status |= gr_mul_two(b, b, ctx);
            status |= gr_mul(a, t, u, ctx);
        }
    }

    GR_TMP_CLEAR2(t, u, ctx);

    return status;
}

int
gr_generic_chebyshev_u_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr a, b;
    fmpz_t n1;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(n))
        return gr_one(y, ctx);

    if (fmpz_is_one(n))
        return gr_mul_two(y, x, ctx);

    if (fmpz_sgn(n) < 0)
    {
        if (fmpz_equal_si(n, -1))
            return gr_zero(y, ctx);

        fmpz_init(n1);
        fmpz_add_ui(n1, n, 2);
        fmpz_neg(n1, n1);
        status = gr_generic_chebyshev_u_fmpz(y, n1, x, ctx);
        status |= gr_neg(y, y, ctx);
        fmpz_clear(n1);
        return status;
    }

    if (gr_is_zero(x, ctx) == T_TRUE)
    {
        int c = fmpz_fdiv_ui(n, 4);
        return gr_set_si(y, (c == 0) - (c == 2), ctx);
    }

    if (gr_is_one(x, ctx) == T_TRUE)
    {
        fmpz_init(n1);
        fmpz_add_ui(n1, n, 1);
        status |= gr_set_fmpz(y, n1, ctx);
        fmpz_clear(n1);
        return status;
    }

    if (gr_is_neg_one(x, ctx) == T_TRUE)
    {
        fmpz_init(n1);
        fmpz_add_ui(n1, n, 1);
        if (fmpz_is_odd(n))
            fmpz_neg(n1, n1);
        status |= gr_set_fmpz(y, n1, ctx);
        fmpz_clear(n1);
        return status;
    }

    GR_TMP_INIT2(a, b, ctx);
    fmpz_init(n1);
    fmpz_tdiv_q_2exp(n1, n, 1);

    status |= gr_generic_chebyshev_u2_fmpz(a, b, n1, x, ctx);

    if (fmpz_is_even(n))
    {
        status |= gr_add(y, a, b, ctx);
        status |= gr_sub(b, a, b, ctx);
        status |= gr_mul(y, y, b, ctx);
    }
    else
    {
        status |= gr_submul(b, a, x, ctx);
        status |= gr_mul(y, a, b, ctx);
        status |= gr_mul_two(y, y, ctx);
        status |= gr_neg(y, y, ctx);
    }

    GR_TMP_CLEAR2(a, b, ctx);
    fmpz_clear(n1);
    return status;
}

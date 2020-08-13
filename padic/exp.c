/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

slong _padic_exp_bound(slong v, slong N, const fmpz_t p)
{
    if (fmpz_fits_si(p))
    {
        fmpz_t n, d, f;
        slong i;

        fmpz_init(n);
        fmpz_init(d);
        fmpz_init(f);

        fmpz_sub_ui(f, p, 1);
        fmpz_mul_ui(n, f, N);
        fmpz_sub_ui(n, n, 1);
        fmpz_mul_ui(d, f, v);
        fmpz_sub_ui(d, d, 1); 
        fmpz_cdiv_q(f, n, d);
        i = fmpz_get_si(f);

        fmpz_clear(n);
        fmpz_clear(d);
        fmpz_clear(f);

        return i;
    }
    else
    {
        return (N + (v - 1)) / v;
    }
}

void _padic_exp(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N)
{
    if (N < 1024)
        _padic_exp_rectangular(rop, u, v, p, N);
    else
        _padic_exp_balanced(rop, u, v, p, N);
}

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const slong N  = padic_prec(rop);
    const slong v  = padic_val(op);
    const fmpz *p = ctx->p;

    if (padic_is_zero(op))
    {
        padic_one(rop);
        return 1;
    }

    if ((fmpz_equal_ui(p, 2) && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
        {
            _padic_exp(padic_unit(rop), padic_unit(op), padic_val(op), p, N);
            padic_val(rop) = 0;
        }
        else
        {
            padic_one(rop);
        }
        return 1;
    }
}


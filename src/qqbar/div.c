/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "qqbar.h"

void
qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(y))
    {
        flint_throw(FLINT_ERROR, "qqbar_div: division by zero\n");
    }
    else if (qqbar_is_zero(x))
    {
        qqbar_zero(res);
    }
    else if (qqbar_is_one(x))
    {
        qqbar_inv(res, y);
    }
    else if (qqbar_is_one(y))
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_neg_one(x))
    {
        qqbar_inv(res, y);
        qqbar_neg(res, res);
    }
    else if (qqbar_is_neg_one(y))
    {
        qqbar_neg(res, x);
    }
    else if (qqbar_is_rational(y))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        _qqbar_get_fmpq(c, a, y);
        qqbar_scalar_op(res, x, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (qqbar_is_rational(x))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        _qqbar_get_fmpq(a, c, x);
        qqbar_inv(res, y);
        qqbar_scalar_op(res, res, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (_qqbar_fast_detect_simple_principal_surd(x) &&
             _qqbar_fast_detect_simple_principal_surd(y))
    {
        /* (p/q)^(1/d) * (r/s)^(1/e) */

        fmpq_t t, u;
        ulong d, e, f, g;

        d = qqbar_degree(x);
        e = qqbar_degree(y);
        g = n_gcd(d, e);
        f = (d / g) * e;

        fmpq_init(t);
        fmpq_init(u);

        fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(x));
        fmpz_set(fmpq_denref(t), QQBAR_COEFFS(x) + d);
        fmpz_neg(fmpq_numref(u), QQBAR_COEFFS(y));
        fmpz_set(fmpq_denref(u), QQBAR_COEFFS(y) + e);

        fmpq_pow_si(t, t, e / g);
        fmpq_pow_si(u, u, d / g);
        fmpq_div(t, t, u);

        /* todo: recycle the existing enclosures instead of computing
           a numerical f-th root from scratch */
        qqbar_fmpq_root_ui(res, t, f);

        fmpq_clear(t);
        fmpq_clear(u);
    }
    else
    {
        qqbar_binary_op(res, x, y, 3);
    }
}

void
qqbar_div_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, y);
    qqbar_div(res, x, t);
    qqbar_clear(t);
}

void
qqbar_div_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, y);
    qqbar_div(res, x, t);
    qqbar_clear(t);
}

void
qqbar_div_ui(qqbar_t res, const qqbar_t x, ulong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, y);
    qqbar_div(res, x, t);
    qqbar_clear(t);
}

void
qqbar_div_si(qqbar_t res, const qqbar_t x, slong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, y);
    qqbar_div(res, x, t);
    qqbar_clear(t);
}

void
qqbar_fmpq_div(qqbar_t res, const fmpq_t x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, x);
    qqbar_div(res, t, y);
    qqbar_clear(t);
}

void
qqbar_fmpz_div(qqbar_t res, const fmpz_t x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, x);
    qqbar_div(res, t, y);
    qqbar_clear(t);
}

void
qqbar_ui_div(qqbar_t res, ulong x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, x);
    qqbar_div(res, t, y);
    qqbar_clear(t);
}

void
qqbar_si_div(qqbar_t res, slong x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, x);
    qqbar_div(res, t, y);
    qqbar_clear(t);
}


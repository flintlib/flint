/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_div(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)
{
    if (ca_qqbar_is_zero(y))
    {
        flint_printf("ca_qqbar_div: division by zero\n");
        flint_abort();
    }
    else if (ca_qqbar_is_zero(x))
    {
        ca_qqbar_zero(res);
    }
    else if (ca_qqbar_is_one(x))
    {
        ca_qqbar_inv(res, y);
    }
    else if (ca_qqbar_is_one(y))
    {
        ca_qqbar_set(res, x);
    }
    else if (ca_qqbar_is_neg_one(x))
    {
        ca_qqbar_inv(res, y);
        ca_qqbar_neg(res, res);
    }
    else if (ca_qqbar_is_neg_one(y))
    {
        ca_qqbar_neg(res, x);
    }
    else if (ca_qqbar_is_rational(y))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(c, CA_QQBAR_COEFFS(y));
        fmpz_set(a, CA_QQBAR_COEFFS(y) + 1);

        ca_qqbar_scalar_op(res, x, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else if (ca_qqbar_is_rational(x))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_neg(a, CA_QQBAR_COEFFS(x));
        fmpz_set(c, CA_QQBAR_COEFFS(x) + 1);

        ca_qqbar_inv(res, y);
        ca_qqbar_scalar_op(res, res, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else
    {
        ca_qqbar_binary_op(res, x, y, 3);
    }
}

void
ca_qqbar_div_fmpq(ca_qqbar_t res, const ca_qqbar_t x, const fmpq_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_fmpq(t, y);
    ca_qqbar_div(res, x, t);
    ca_qqbar_clear(t);
}

void
ca_qqbar_div_fmpz(ca_qqbar_t res, const ca_qqbar_t x, const fmpz_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_fmpz(t, y);
    ca_qqbar_div(res, x, t);
    ca_qqbar_clear(t);
}

void
ca_qqbar_div_ui(ca_qqbar_t res, const ca_qqbar_t x, ulong y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_ui(t, y);
    ca_qqbar_div(res, x, t);
    ca_qqbar_clear(t);
}

void
ca_qqbar_div_si(ca_qqbar_t res, const ca_qqbar_t x, slong y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_si(t, y);
    ca_qqbar_div(res, x, t);
    ca_qqbar_clear(t);
}

void
ca_qqbar_fmpq_div(ca_qqbar_t res, const fmpq_t x, const ca_qqbar_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_fmpq(t, x);
    ca_qqbar_div(res, t, y);
    ca_qqbar_clear(t);
}

void
ca_qqbar_fmpz_div(ca_qqbar_t res, const fmpz_t x, const ca_qqbar_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_fmpz(t, x);
    ca_qqbar_div(res, t, y);
    ca_qqbar_clear(t);
}

void
ca_qqbar_ui_div(ca_qqbar_t res, ulong x, const ca_qqbar_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_ui(t, x);
    ca_qqbar_div(res, t, y);
    ca_qqbar_clear(t);
}

void
ca_qqbar_si_div(ca_qqbar_t res, slong x, const ca_qqbar_t y)
{
    ca_qqbar_t t;
    ca_qqbar_init(t);
    ca_qqbar_set_si(t, x);
    ca_qqbar_div(res, t, y);
    ca_qqbar_clear(t);
}


/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(x))
    {
        qqbar_neg(res, y);
    }
    else if (qqbar_is_zero(y))
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_rational(y))
    {
        fmpz_t a, b, c;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);

        fmpz_set(a, QQBAR_COEFFS(y) + 1);
        fmpz_set(b, QQBAR_COEFFS(y));
        fmpz_set(c, QQBAR_COEFFS(y) + 1);

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

        fmpz_neg(a, QQBAR_COEFFS(x) + 1);
        fmpz_neg(b, QQBAR_COEFFS(x));
        fmpz_set(c, QQBAR_COEFFS(x) + 1);

        qqbar_scalar_op(res, y, a, b, c);

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
    }
    else
    {
        qqbar_binary_op(res, x, y, 1);
    }
}

void
qqbar_sub_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, y);
    qqbar_sub(res, x, t);
    qqbar_clear(t);
}

void
qqbar_sub_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, y);
    qqbar_sub(res, x, t);
    qqbar_clear(t);
}

void
qqbar_sub_ui(qqbar_t res, const qqbar_t x, ulong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, y);
    qqbar_sub(res, x, t);
    qqbar_clear(t);
}

void
qqbar_sub_si(qqbar_t res, const qqbar_t x, slong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, y);
    qqbar_sub(res, x, t);
    qqbar_clear(t);
}

void
qqbar_fmpq_sub(qqbar_t res, const fmpq_t x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, x);
    qqbar_sub(res, t, y);
    qqbar_clear(t);
}

void
qqbar_fmpz_sub(qqbar_t res, const fmpz_t x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, x);
    qqbar_sub(res, t, y);
    qqbar_clear(t);
}

void
qqbar_ui_sub(qqbar_t res, ulong x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, x);
    qqbar_sub(res, t, y);
    qqbar_clear(t);
}

void
qqbar_si_sub(qqbar_t res, slong x, const qqbar_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, x);
    qqbar_sub(res, t, y);
    qqbar_clear(t);
}

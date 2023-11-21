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
qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y)
{
    if (qqbar_is_zero(x))
    {
        qqbar_set(res, y);
    }
    else if (qqbar_is_zero(y))
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_rational(y))
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_init(b);

        _qqbar_get_fmpq(b, a, y);
        qqbar_scalar_op(res, x, a, b, a);

        fmpz_clear(a);
        fmpz_clear(b);
    }
    else if (qqbar_is_rational(x))
    {
        fmpz_t a, b;

        fmpz_init(a);
        fmpz_init(b);

        _qqbar_get_fmpq(b, a, x);
        qqbar_scalar_op(res, y, a, b, a);

        fmpz_clear(a);
        fmpz_clear(b);
    }
    else
    {
        qqbar_binary_op(res, x, y, 0);
    }
}

void
qqbar_add_fmpq(qqbar_t res, const qqbar_t x, const fmpq_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpq(t, y);
    qqbar_add(res, x, t);
    qqbar_clear(t);
}

void
qqbar_add_fmpz(qqbar_t res, const qqbar_t x, const fmpz_t y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_fmpz(t, y);
    qqbar_add(res, x, t);
    qqbar_clear(t);
}

void
qqbar_add_ui(qqbar_t res, const qqbar_t x, ulong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_ui(t, y);
    qqbar_add(res, x, t);
    qqbar_clear(t);
}

void
qqbar_add_si(qqbar_t res, const qqbar_t x, slong y)
{
    qqbar_t t;
    qqbar_init(t);
    qqbar_set_si(t, y);
    qqbar_add(res, x, t);
    qqbar_clear(t);
}


/*
    Copyright (C) 2012-2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_add_error_arf(arb_t x, const arf_t err)
{
    mag_t t;

    if (arf_is_zero(err))
        return;

    if (mag_is_zero(arb_radref(x)))
    {
        arf_get_mag(arb_radref(x), err);
        return;
    }

    mag_init(t);
    arf_get_mag(t, err);
    mag_add(arb_radref(x), arb_radref(x), t);
    mag_clear(t);
}

void
arb_add_error_2exp_si(arb_t x, slong err)
{
    fmpz_t t;

    if (mag_is_zero(arb_radref(x)))
    {
        mag_one(arb_radref(x));
        mag_mul_2exp_si(arb_radref(x), arb_radref(x), err);
        return;
    }

    fmpz_init(t);
    fmpz_set_si(t, err);
    mag_add_2exp_fmpz(arb_radref(x), arb_radref(x), t);
    fmpz_clear(t);
}

void
arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err)
{
    if (mag_is_zero(arb_radref(x)))
    {
        mag_one(arb_radref(x));
        mag_mul_2exp_fmpz(arb_radref(x), arb_radref(x), err);
        return;
    }

    mag_add_2exp_fmpz(arb_radref(x), arb_radref(x), err);
}

void
arb_add_error(arb_t x, const arb_t err)
{
    mag_t u;

    if (arb_is_zero(err))
        return;

    if (mag_is_zero(arb_radref(x)))
    {
        arb_get_mag(arb_radref(x), err);
        return;
    }

    mag_init(u);
    arb_get_mag(u, err);
    mag_add(arb_radref(x), arb_radref(x), u);
    mag_clear(u);
}


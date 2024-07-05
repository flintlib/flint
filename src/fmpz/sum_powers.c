/*
    Copyright (C) 2024 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void
fmpz_sum_powers_horner(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz_t t;
    fmpz_init_set_ui(t, 1);

    for (ulong i = 1; i < exp; i++)
    {
        fmpz_mul(t, t, g);
        fmpz_add_ui(t, t, 1);
    }

    fmpz_swap(f, t);
    fmpz_clear(t);
}

void
fmpz_sum_powers_div(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz_t t;
    fmpz_init(t);

    // t = g^(e + 1) - 1
    fmpz_pow_ui(t, g, exp + 1);
    fmpz_sub_ui(t, t, 1);

    // f = g - 1
    fmpz_sub_ui(f, g, 1);

    // f = (g^(e + 1) - 1) / (g - 1)
    fmpz_tdiv(f, t, g);

    fmpz_clear(t);
}

void
fmpz_sum_powers(fmpz_t f, const fmpz_t g, ulong exp)
{
    if (exp == 0)
    {
        fmpz_zero(f);
    }
    else if (exp == 1)
    {
        fmpz_add_ui(f, g, 1);
    }
    else if (exp <= 100)
    {
        _fmpz_sum_powers_horner(f, g, exp);
    }
    else
    {
        _fmpz_sum_powers_div(f, g, exp);
    }
}

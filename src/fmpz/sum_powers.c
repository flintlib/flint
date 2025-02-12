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

    if (exp == 0 || fmpz_is_zero(g))
    {
        fmpz_one(f);
    }
    else if (fmpz_is_one(g))
    {
        fmpz_set_ui(f, exp + 1);
    }
    else
    {
        fmpz_init_set_ui(t, 1);

        for (ulong i = 1; i <= exp; i++)
        {
            fmpz_mul(t, t, g);
            fmpz_add_ui(t, t, 1);
        }

        fmpz_swap(f, t);
        fmpz_clear(t);
    }
}

void
fmpz_sum_powers_div(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz_t t;

    if (exp == 0 || fmpz_is_zero(g))
    {
        fmpz_one(f);
    }
    else if (fmpz_is_one(g))
    {
        fmpz_set_ui(f, exp + 1);
    }
    else
    {
        fmpz_init(t);

        /* t = g - 1 */
        fmpz_sub_ui(t, g, 1);

        /* f = g^(e + 1) - 1 */
        fmpz_pow_ui(f, g, exp + 1);
        fmpz_sub_ui(f, f, 1);

        /* f = (g^(e + 1) - 1) / (g - 1) */
        fmpz_tdiv_q(f, f, t);

        fmpz_clear(t);
    }
}

void
fmpz_sum_powers(fmpz_t f, const fmpz_t g, ulong exp)
{
    if (exp == 0 || fmpz_is_zero(g))
    {
        fmpz_one(f);
    }
    else if (exp == 1)
    {
        fmpz_add_ui(f, g, 1);
    }
    else if (exp <= 100)
    {
        fmpz_sum_powers_horner(f, g, exp);
    }
    else
    {
        fmpz_sum_powers_div(f, g, exp);
    }
}

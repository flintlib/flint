/*
    Copyright (C) 2024 Matthias Gessinger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_factor.h"

void
fmpz_sum_divisors(fmpz_t f, const fmpz_t g)
{
    fmpz_factor_t factors;
    fmpz * temp;
    ulong exp;
    slong i;

    fmpz * p;

    if (fmpz_is_zero(g) || fmpz_is_pm1(g))
    {
        fmpz_set(f, g);
    }
    else
    {
        fmpz_factor_init(factors);

        fmpz_factor(factors, g);
        i = 0;

        temp = factors->p + i;
        exp = factors->exp[i];

        fmpz_sum_powers(temp, temp, exp);

        for (i = 1; i < factors->num; i++)
        {
            p = factors->p + i;
            exp = factors->exp[i];

            fmpz_sum_powers(p, p, exp);

            fmpz_mul(temp, temp, p);
        }

        if (factors->sign == -1)
        {
            fmpz_neg(f, temp);
        }
        else
        {
            fmpz_swap(f, temp);
        }

        fmpz_factor_clear(factors);
    }
}

void
fmpz_sum_divisors_proper(fmpz_t f, const fmpz_t g)
{
    fmpz_t temp;
    fmpz_init(temp);

    fmpz_sum_divisors(temp, g);

    fmpz_sub(f, temp, g);

    fmpz_clear(temp);
}

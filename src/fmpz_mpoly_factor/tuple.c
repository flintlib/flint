/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"


void tuple_print(fmpz * alpha, slong n)
{
    slong j;
    for (j = 0; j < n; j++)
    {
        fmpz_print(alpha + j);
        flint_printf(j + 1 < n ? ", " : "\n");
    }
}

/* ensure that the first m values change upon the next call to tuple_next*/
void tuple_saturate(fmpz * alpha, slong n, slong m)
{
    slong i;

    for (i = m + 1; i < n; i++)
    {
        fmpz_add(alpha + m, alpha + m, alpha + i);
        fmpz_zero(alpha + i);
    }

    if (m < n && fmpz_is_zero(alpha + m))
    {
        for (i = 0; i < m; i++)
            if (!fmpz_is_zero(alpha + i))
                return;
        fmpz_one(alpha + m);
    }
}


void tuple_next(fmpz * alpha, slong n)
{
    slong i, t1, t2, t3;
    fmpz_t sum;

    fmpz_init(sum);
    for (i = 0; i < n; i++)
        fmpz_add(sum, sum, alpha + i);

    i = n - 1;
    while(i >= 0 && fmpz_is_zero(alpha + i))
        i--;
    t1 = i;
    while(i >= 0 && fmpz_cmp(alpha + i, sum) != 0)
        i--;
    t2 = i;
    while(i >= 0 && fmpz_cmp(alpha + i, sum) == 0)
        i--;
    t3 = i;

    if (t1 > 0 && t1 != t2)
    {
        fmpz_swap(alpha + t1, alpha + n - 1);
        fmpz_sub_ui(alpha + n - 1, alpha + n - 1, 1);
        fmpz_add_ui(alpha + t1 - 1, alpha + t1 - 1, 1);
    }
    else if (t1 > 0 && t1 == t2 && t3 >= 0)
    {
        fmpz_add_ui(alpha + t3, alpha + t3, 1);
        fmpz_zero(alpha + t3 + 1);
        fmpz_sub_ui(alpha + n - 1, sum, 1);
    }
    else
    {
        fmpz_add_ui(alpha + n - 1, alpha + 0, 1);
        if (n > 1)
            fmpz_zero(alpha + 0);
    }

    fmpz_clear(sum);
}


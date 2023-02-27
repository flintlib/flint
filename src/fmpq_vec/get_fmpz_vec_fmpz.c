/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_vec.h"

void
_fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den, const fmpq * a, slong len)
{
    slong i;

    if (len < 1)
    {
        fmpz_one(den);
        return;
    }
    else if (len == 1)
    {
        fmpz_set(num + 0, fmpq_numref(a + 0));
        fmpz_set(den, fmpq_denref(a + 0));
        return;
    }

    fmpz_lcm(den, fmpq_denref(a + 0), fmpq_denref(a + 1));
    for (i = 2; i < len; i++)
        fmpz_lcm(den, den, fmpq_denref(a + i));

    if (fmpz_is_one(den))
    {
        for (i = 0; i < len; i++)
            fmpz_set(num + i, fmpq_numref(a + i));
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            fmpz_divexact(num + i, den, fmpq_denref(a + i));
            fmpz_mul(num + i, num + i, fmpq_numref(a + i));
        }
    }
}

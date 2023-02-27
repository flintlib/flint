/*
    Copyright (C) 2015, Elias Tsigaridas
    Copyright (C) 2016, Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void _fmpz_poly_scale_2exp(fmpz * pol, slong len, slong k)
{
    slong i;
    ulong p, z;

    if (k == 0 || len <= 1)
        return;
    else if (k > 0)
    {
        if (fmpz_is_zero(pol))
            z = ULONG_MAX;
        else
            z = fmpz_val2(pol);

        p = k;
        for (i = 1; i < len; i++, p += k)
        {
            if (!fmpz_is_zero(pol + i))
            {
                z = FLINT_MIN(z, fmpz_val2(pol + i) + p);
                fmpz_mul_2exp(pol + i, pol + i, p);
            }
        }
    }
    else
    {
        if (fmpz_is_zero(pol + len - 1))
            z = ULONG_MAX;
        else
            z = fmpz_val2(pol + len - 1);

        p = -k;
        for (i = len - 2; i >= 0; i--, p -= k)
        {
            if (!fmpz_is_zero(pol + i))
            {
                z = FLINT_MIN(z, fmpz_val2(pol + i) + p);
                fmpz_mul_2exp(pol + i, pol + i, p);
            }
        }
    }

    if (z)
    {
        for (i = 0; i < len; i++)
            fmpz_fdiv_q_2exp(pol + i, pol + i, z);
    }
}

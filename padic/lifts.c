/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

slong * _padic_lifts_exps(slong *n, slong N)
{
    slong *a, i;

    *n = FLINT_CLOG2(N) + 1;

    a = flint_malloc((*n) * sizeof(slong));
    for (a[i = 0] = N; a[i] > 1; i++)
        a[i + 1] = (a[i] + 1) / 2;

    return a;
}

void _padic_lifts_pows(fmpz *pow, const slong *a, slong n, const fmpz_t p)
{
    if (n == 1)
    {
        fmpz_set(pow + 0, p);
    }
    else  /* n > 1 */
    {
        slong i = n - 1;
        fmpz_t t = {WORD(1)};

        fmpz_set(pow + i, p);

        for (i--; i >= 1; i--)
        {
            if (a[i] & WORD(1))
            {
                fmpz_mul(pow + i, t, pow + (i + 1));
                fmpz_mul(t, t, t);
            }
            else
            {
                fmpz_mul(t, t, pow + (i + 1));
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            }
        }
        {
            if (a[i] & WORD(1))
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }
        fmpz_clear(t);
    }
}


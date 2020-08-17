/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arith.h"

void
arith_bell_number_vec_recursive(fmpz * b, slong n)
{
    slong i, k;
    fmpz * t;

    if (n < BELL_NUMBER_TAB_SIZE)
    {
        for (i = 0; i < n; i++)
            fmpz_set_ui(b + i, bell_number_tab[i]);
        return;
    }

    n -= 1;
    t = _fmpz_vec_init(n);

    fmpz_one(t);
    fmpz_one(b);
    fmpz_one(b + 1);

    for (i = 1; i < n; i++)
    {
        fmpz_set(t + i, t);
        for (k = i; k > 0; k--)
            fmpz_add(t + k - 1, t + k - 1, t + k);
        fmpz_set(b + i + 1, t);
    }

    _fmpz_vec_clear(t, n);
}

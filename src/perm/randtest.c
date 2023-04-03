/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "flint.h"
#include "ulong_extras.h"

int _perm_randtest(slong *vec, slong n, flint_rand_t state)
{
    slong i, j, t;

    int parity = 0;

    for (i = 0; i < n; i++)
        vec[i] = i;

    for (i = n - 1; i > 0; i--)
    {
        j = n_randint(state, i + 1);
        parity ^= (i != j);
        t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }

    return parity;
}

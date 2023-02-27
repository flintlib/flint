/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"

int
_perm_parity(const slong *vec, slong n)
{
    slong i, k;
    int * encountered;
    int parity;
    TMP_INIT;

    if (n <= 1)
        return 0;

    TMP_START;

    parity = 0;
    encountered = (int *) TMP_ALLOC(n*sizeof(int));

    for (i = 0; i < n; i++)
       encountered[i] = 0;

    for (i = 0; i < n; i++)
    {
        if (encountered[i] != 0)
        {
            parity ^= 1;
        }
        else
        {
            k = i;
            do
            {
                k = vec[k];
                encountered[k] = 1;
            }
            while (k != i);
        }
    }

    TMP_END;

    return parity;
}

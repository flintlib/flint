/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

/* assume mod is small so no overflow */
ulong
dlog_table_init(dlog_table_t t, ulong a, ulong mod)
{
    int k;
    ulong ak;
    t->mod = mod;
    t->table = flint_malloc(mod * sizeof(ulong));
    ak = 1; k = 0;

    /* warning: do not check a is invertible modulo mod */
    do
    {
        t->table[ak] = k++;
        ak = (ak * a) % mod;
    }

    while (ak != 1);
    return 1;
}

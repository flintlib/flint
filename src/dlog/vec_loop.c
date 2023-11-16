/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "dlog.h"

/* vector of log(k,a)*loga % order in Z/modZ */
void
dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong x, vx;
    dlog_vec_fill(v, nv, DLOG_NOT_FOUND);
    x = 1; vx = 0;
    do
    {
        if (x < nv)
            v[x] = vx;
        x = nmod_mul(x, a, mod);
        vx = nmod_add(vx, va, order);
    }
    while (x != 1);
    for (x = mod.n + 1; x < nv; x++)
        v[x] = v[x - mod.n];
}

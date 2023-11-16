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
dlog_vec_loop_add(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong x, xp, vx;
    vx = 0;
    for (x = a; x != 1; x = nmod_mul(x, a, mod))
    {
        vx = nmod_add(vx, va, order);
        for(xp = x; xp < nv; xp+=mod.n)
            if (v[xp] != DLOG_NONE)
                v[xp] = nmod_add(v[xp], vx, order);
    }
}

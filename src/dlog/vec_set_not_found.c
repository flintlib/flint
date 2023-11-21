/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

void
dlog_vec_set_not_found(ulong *v, ulong nv, nmod_t mod)
{
    n_factor_t fac;
    ulong i;

    n_factor_init(&fac);
    n_factor(&fac, mod.n, 1);
    for (i = 0; i < fac.num; i++)
    {
        ulong p, k;
        p = fac.p[i];
        for (k = p; k < nv; k += p)
            v[k] = DLOG_NOT_FOUND;
    }
}

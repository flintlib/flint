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

void
dlog_vec_sieve_add(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong * w, k;
    /* store size */
    w = flint_malloc(nv * sizeof(ulong));
    dlog_vec_sieve(w, nv, a, va, mod, na, order);
    /* write in v */
    for (k = 0; k < nv; k++)
        if (v[k] != DLOG_NOT_FOUND)
            v[k] = nmod_add(v[k], w[k], order);
    flint_free(w);
}

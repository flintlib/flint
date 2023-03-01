/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dlog.h"

void
dlog_vec_sieve_add_precomp(ulong *v, ulong nv, dlog_precomp_t pre, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong * w, k;
    /* store size */
    w = flint_malloc(nv * sizeof(ulong));
    dlog_vec_sieve_precomp(w, nv, pre, a, va, mod, na, order);
    /* write in v */
    for (k = 0; k < nv; k++)
        if (v[k] != DLOG_NOT_FOUND)
            v[k] = nmod_add(v[k], w[k], order);
    flint_free(w);
}

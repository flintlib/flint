/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod.h"
#include "dlog.h"

int
apow_cmp(const apow_t * x, const apow_t * y)
{
    return (x->ak < y->ak) ? -1 : (x->ak > y->ak);
}

ulong
dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong mod, ulong n, ulong m)
{
    ulong k, ak;
    if (m >= n) m = n;
    t->table = (apow_t *)flint_malloc(m * sizeof(apow_t));

    nmod_init(&t->mod, mod);
    t->m = m;
    t->g = n / m + 1;

    for (k = 0, ak = 1; k < m; k++)
    {
        t->table[k].k = k;
        t->table[k].ak = ak;
        ak = nmod_mul(ak, a, t->mod);
    }

    t->am = nmod_inv(ak, t->mod);
    qsort(t->table, m, sizeof(apow_t), (int(*)(const void*,const void*))apow_cmp);
    return t->g;
}

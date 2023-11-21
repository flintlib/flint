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

ulong
dlog_power_init(dlog_power_t t, ulong a, ulong mod, ulong p, ulong e, ulong num)
{
    int k;
    nmod_init(&t->mod, mod);
    t->p = p;
    t->e = e;

    t->apk = flint_malloc(e * sizeof(ulong));
    t->pre = flint_malloc(sizeof(dlog_precomp_struct));

    t->apk[0] = nmod_inv(a, t->mod);
    for (k = 1; k < e; k++)
        t->apk[k] = nmod_pow_ui(t->apk[k-1], p, t->mod);

    dlog_precomp_p_init(t->pre, nmod_inv(t->apk[e-1], t->mod), mod, p, e * num);
    return e * t->pre->cost;
}

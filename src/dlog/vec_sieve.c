/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

#define vbs 0

/* TODO: tune the limit dlog -> index calculus */
void
dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong p1 = 50; /* FIXME: tune this limit! */
    dlog_precomp_t pre;

    dlog_precomp_n_init(pre, a, mod.n, na, p1);
    dlog_vec_sieve_precomp(v, nv, pre, a, va, mod, na, order);
    dlog_precomp_clear(pre);
}

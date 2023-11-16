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
dlog_vec(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    if (va == 0)
        return;
    if (na * DLOG_LOOP_MAX_FACTOR < nv)
        dlog_vec_loop(v, nv, a, va, mod, na, order);
    else
        dlog_vec_sieve(v, nv, a, va, mod, na, order);
}

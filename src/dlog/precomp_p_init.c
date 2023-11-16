/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

/* we known the order is prime */
void
dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num)
{
    if (p < DLOG_TABLE_P_LIM)
    {
        dlog_precomp_small_init(pre, a, mod, p, num);
    }
    else
    {
        ulong m;
        m = dlog_bsgs_size(p, num);
        pre->type = DLOG_BSGS;
        pre->cost = dlog_bsgs_init(pre->t.bsgs, a, mod, p, m);
    }
}

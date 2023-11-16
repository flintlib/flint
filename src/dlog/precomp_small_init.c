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
dlog_precomp_small_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)
{
    if (n <= 3)
    {
        pre->type = DLOG_23;
        pre->cost = dlog_order23_init(pre->t.order23, a);
    }
    else
    {
        if (mod < DLOG_TABLE_LIM)
        {
            pre->type = DLOG_TABLE;
            pre->cost = dlog_table_init(pre->t.table, a, mod);
        }
        else
        {
            pre->type = DLOG_BSGS;
            pre->cost = dlog_bsgs_init(pre->t.bsgs, a, mod, n, n);
        }
    }
}

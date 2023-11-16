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
dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num)
{
    if (pe < DLOG_TABLE_PE_LIM)
    {
        dlog_precomp_small_init(pre, a, mod, pe, num);
    }
    else
    {
        if (e == 1)
        {
            dlog_precomp_p_init(pre, a, mod, p, num);
        }
        else
        {
            pre->type = DLOG_POWER;
            pre->cost = dlog_power_init(pre->t.power, a, mod, p, e, num);
        }
    }
}

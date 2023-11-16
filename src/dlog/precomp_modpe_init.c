/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

/* log mod p^e */
void
dlog_precomp_modpe_init(dlog_precomp_t pre, ulong a, ulong p, ulong e, ulong pe, ulong num)
{
    if (pe < DLOG_TABLE_MODPE_LIM)
    {
        dlog_precomp_small_init(pre, a, pe, pe - pe / p, num);
        return;
    }
    else
    {
        if (e > 1)
        {
            pre->type = DLOG_MODPE;
            pre->cost = dlog_modpe_init(pre->t.modpe, a, p, e, pe, num);
        }
        else
        {
            dlog_precomp_n_init(pre, a, p, p - 1, num);
        }
    }
}

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
dlog_1modpe(const dlog_1modpe_t t, ulong b1, ulong p, ulong e, nmod_t pe)
{
    if (e == 1)
        return 0;
    else
    {
        ulong logb1;
        logb1 = dlog_1modpe_1modp(b1, p, e, t->inv1p, pe);
        /* only need mod p^(e-1) */
        return nmod_mul(logb1, t->invloga1, pe);
    }
}

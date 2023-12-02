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

void
dlog_1modpe_init(dlog_1modpe_t t, ulong a1, ulong p, ulong e, nmod_t pe)
{
    if (e == 1)
    {
        t->inv1p = 1;
        t->invloga1 = 0;
    }
    else
    {
        ulong loga1;
        if (a1 == 1)
            flint_throw(FLINT_ERROR, __func__);
        t->inv1p = nmod_inv(1 + p, pe); /* 1 - p + p^2 - ... */
        loga1 = dlog_1modpe_1modp(a1, p, e, t->inv1p, pe);
        /* only need inverse mod p^(e-1) but does not hurt */
        t->invloga1 = nmod_inv(loga1, pe);
    }
}

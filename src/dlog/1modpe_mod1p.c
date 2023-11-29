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

/* for odd prime p, assume b1 = 1 mod p */
ulong
dlog_1modpe_1modp(ulong b1, ulong p, ulong e, ulong inv1p, nmod_t pe)
{
    int f;
    ulong x, xf, pf, pf1;
    pf1 = 1;
    pf = p;
    x = 0;
    for (f = 1; f < e; f++)
    {
        if (b1 % pf != 1)
        {
            flint_throw(FLINT_ERROR, "ERROR dlog_1modpe_1modp: %wu %% %wu != 1 mod %wu\n\n", b1, pf, pe.n);
        }
        xf = (b1 - 1) / pf;
        xf = (xf % p) * pf1;
        x += xf;
        b1 = nmod_mul(b1, nmod_pow_ui(inv1p, xf, pe), pe);
        pf1 = pf;
        pf *= p;
    }
    return x;
}

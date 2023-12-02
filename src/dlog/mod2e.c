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
dlog_mod2e_1mod4(ulong b1, ulong e, ulong inv5, nmod_t pe)
{
    slong f;
    ulong pf1, pf, x, xf;
    pf1 = 1;
    pf = 4;
    x = 0;
    for (f = 2; f < e; f++)
    {
        if (b1 % pf != 1)
        {
            flint_throw(FLINT_ERROR, "ERROR dlog_mod2e: %wu %% %wu != 1 mod %wu\n\n", b1, pf, pe.n);
        }
        xf = (b1 - 1) / pf;
        xf = (f == 2) ? xf % 4 : (xf % 2) * (pf1 / 2);
        b1 = nmod_mul(b1, nmod_pow_ui(inv5, xf, pe), pe);
        x += xf;
        pf1 = pf;
        pf *= 2;
    }
    return x;
}

ulong
dlog_mod2e(const dlog_modpe_t t, ulong b1)
{
    if (t->e == 2)
        return (b1 % 4) == 3;
    else
        return dlog_mod2e_1mod4(b1, t->e, t->inva, t->pe);
}

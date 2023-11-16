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
dlog_modpe(const dlog_modpe_t t, ulong b)
{
    ulong x;
    if (t->p == 2)
        return dlog_mod2e(t, b);

    x = dlog_precomp(t->modp, b % t->p);

    if (t->e > 1)
    {
        ulong b1, y;
#if 0
        b1 = nmod_mul(b, nmod_pow_ui(t->inva, x, t->pe), t->pe);
        y = dlog_1modpe(t->modpe.rec, b1, t->p, t->e, t->pe);
        y = y % t->pe1;
        x = x + (t->p - 1) * y;
#else
        b1 = nmod_pow_ui(b, t->p - 1, t->pe);
        y = dlog_1modpe(t->modpe, b1, t->p, t->e, t->pe);
        y = y % t->pe1;
        x = n_submod(x, y % (t->p - 1), t->p - 1);
        x = y + t->pe1 * x;
#endif
    }

    return x;
}

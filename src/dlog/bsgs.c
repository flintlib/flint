/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod.h"
#include "dlog.h"

ulong
dlog_bsgs(const dlog_bsgs_t t, ulong b)
{
    ulong i;
    apow_t c, * x;

    c.ak = b;
    for (i = 0; i < t->g; i++)
    {
        x = bsearch(&c, t->table, t->m, sizeof(apow_t),
            (int(*)(const void*,const void*))apow_cmp);
        if (x != NULL)
            return i * t->m + x->k;
        c.ak = nmod_mul(c.ak, t->am, t->mod);
    }

    flint_throw(FLINT_ERROR, "Exception (dlog_bsgs).  discrete log not found.\n"
            "   table size %wu, cosize %wu mod %wu. %wu not found (a^-m=%wu)\n",
            t->m, t->g, t->mod.n, b, t->am);
}

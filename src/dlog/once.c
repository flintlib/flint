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
dlog_once(ulong b, ulong a, const nmod_t mod, ulong n)
{
    if (b == 1)
        return 0;
    if (n < 50)
    {
        slong k;
        ulong ak = 1;

        for (k = 0; k < n; k++)
        {
            if (ak == b)
                return k;
            ak = nmod_mul(ak, a, mod);
        }

        flint_throw(FLINT_ERROR, "(dlog_once): log(%wu,%wu) mod %wu not found (size %wu)\n", b, a, mod.n, n);
    }
    else
    {
        ulong l;
        dlog_precomp_t pre;
        dlog_precomp_n_init(pre, a, mod.n, n, 1);
        l = dlog_precomp(pre, b);
        dlog_precomp_clear(pre);
        return l;
    }
}

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
dlog_crt(const dlog_crt_t t, ulong b)
{
    int k;
    ulong r = 0;
    for (k = 0; k < t->num; k++)
    {
        ulong bk, rk;
        bk = nmod_pow_ui(b, t->expo[k], t->mod);
        rk = dlog_precomp(t->pre + k, bk);
#if 0
        flint_printf("##[crt-%d]: log(%wu)=log(%wu^%wu) = %wu [size %wu mod %wu]\n",
                k, bk, b, t->expo[k], rk, t->n/t->expo[k], t->mod);
#endif
        r = nmod_add(r, nmod_mul(t->crt_coeffs[k], rk, t->n), t->n);
    }
    return r;
}

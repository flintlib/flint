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
dlog_power(const dlog_power_t t, ulong b)
{
    int k;
    ulong x, pk[30]; /* 3^30*2+1, 2^30*3+1 are primes */

    pk[0] = 1;

    for (k = 1; k < t->e; k++)
       pk[k] = pk[k-1] * t->p;

    x = 0;

    for(k = 0; k < t->e; k++)
    {
      ulong bk, xk;
      bk = nmod_pow_ui(b, pk[t->e-1-k], t->mod);
      xk = dlog_precomp(t->pre, bk);
      b = nmod_mul(b, nmod_pow_ui(t->apk[k], xk, t->mod), t->mod);
      x += xk * pk[k]; /* cannot overflow */
    }

    return x;
}

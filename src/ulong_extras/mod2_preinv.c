/*
    Copyright (C) 2009, 2010, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"

/* Barrett step of NMOD_RED (doc: algorithms/arithmetic/barrett#nmod) */

ulong
n_mod2_preinv(ulong a, ulong n, ulong ninv)
{
    nmod_t mod;
    ulong r;

    FLINT_ASSERT(n != 0);

    mod.n = n;
    mod.ninv = ninv;
    mod.norm = flint_clz(n);
    NMOD_RED(r, a, mod);
    return r;
}

/* NMOD_RED without its final correction. For norm >= 1, the quotient
   estimate before the floor is below a / n by less than 2 / n (proof of
   algorithms/arithmetic/barrett#thm-nmod-red), so the result is below
   n + 2; for norm = 0 it is exact. */
ulong
n_mod2_preinv_lazy(ulong a, ulong n, ulong ninv)
{
    unsigned int norm = flint_clz(n);

    FLINT_ASSERT(n != 0);

    if (norm == 0)
        return (a >= n) ? a - n : a;
    else
    {
        const ulong m = (ninv >> 1) | (UWORD(1) << (FLINT_BITS - 1));
        return a - (n_mulhi(a, m) >> (FLINT_BITS - 1 - norm)) * n;
    }
}

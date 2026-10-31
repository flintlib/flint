/*
    Copyright (C) 2009, 2013, 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod/pow_impl.h"

ulong
n_powmod2_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv)
{
    nmod_t mod;

    FLINT_ASSERT(n != 0);

    if (a >= n)
        a = n_mod2_preinv(a, n, ninv);

    mod.n = n;
    mod.ninv = ninv;
    mod.norm = flint_clz(n);

    return _pow_ui(a, exp, mod);
}

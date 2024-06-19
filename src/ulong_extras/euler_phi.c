/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_euler_phi(ulong n)
{
    int i;
    ulong phi;
    n_factor_t fac;

    if (n < 2)
        return n;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    phi = UWORD(1);
    for (i = 0; i < fac.num; i++)
        phi *= (fac.p[i] - 1) * n_pow(fac.p[i], fac.exp[i] - 1);

    return phi;
}

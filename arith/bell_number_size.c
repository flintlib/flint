/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"

double
arith_bell_number_size(ulong n)
{
    double l, ll, u;

    if (n <= 1)
        return 0;

    /* Using de Bruijn's asymptotic expansion. Not sure if this is an
       upper bound for all n, but suffices at least for n < 2^64. */
    l = log((double) n);
    ll = log(l);
    u = 1.0 / l;

    return 1.4426950408889634074 * n * (l - ll - 1.0 + ll * u
        + 1.0 * u + 0.5 * (ll * u) * (ll * u) + 0.25 * ll * u * u) + 2;
}

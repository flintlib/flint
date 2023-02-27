/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "long_extras.h"
#include "ulong_extras.h"

int z_kronecker(slong a, slong n)
{
    ulong sa = FLINT_SIGN_EXT(a);
    ulong ua = FLINT_ABS(a);
    ulong sn = FLINT_SIGN_EXT(n);
    ulong un = FLINT_ABS(n);
    int en;
    unsigned int r;

    if (a == 0)
        return un == 1;

    if (un == 0)
        return ua == 1;

    count_trailing_zeros(en, un);

    /* make denominator positive */
    r = sa & sn;

    /* make denominator odd */
    un >>= en;
    r ^= (ua ^ (ua>>1)) & (2*en);

    if (en > 0 && !(ua & 1))
        return 0; /* a and n both even */

    /* make numerator positive */
    r ^= (sa & un);

    return _n_jacobi_unsigned(ua, un, r);
}


/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arith.h"

double arith_bernoulli_number_size(ulong n)
{
    double x;

    /* |B_n| < 2 */
    if (n <= 14)
        return 1.0;

    x = 2 + (n + 1) * log(n + 1) * 1.44269504088897;  /* 1/log(2) */
    x -= n * 4.0941911703612822; /* log2(2*pi*e) */

    return x + 2;
}

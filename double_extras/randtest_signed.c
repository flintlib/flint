/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "ulong_extras.h"

double
d_randtest_signed(flint_rand_t state, slong minexp, slong maxexp)
{
    double d, t;
    slong exp, kind;
    d = d_randtest(state);
    exp = minexp + n_randint(state, maxexp - minexp + 1);
    t = ldexp(d, exp);
    kind = n_randint(state, 3);
	if (kind == 2)
        return t;
    else if (kind == 1)
        return -t;
    else
        return 0;
}

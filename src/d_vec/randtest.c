/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "d_vec.h"

void
_d_vec_randtest(double *f, flint_rand_t state, slong len, slong minexp,
                slong maxexp)
{
    slong i;

    for (i = 0; i < len; i++)
        f[i] = d_randtest_signed(state, minexp, maxexp);
}

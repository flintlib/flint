/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fixed.h"

/* log 10 = log 2 + log 5, the two from the three-prime Machin-type set
   (46 atanh(1/31) + 34 atanh(1/49) + 20 atanh(1/161) in the plain
   arctangent form, each term from Zuniga's series where that is
   faster), which gives log 2, log 3 and log 5 from the same three
   series and runs them on the available threads. */
void
fball_const_log10(fball_t res, slong n)
{
    fball_struct v[3];
    slong i;

    for (i = 0; i < 3; i++)
        fball_init(v + i);

    _fixed_log_primes_vec_fball(v, 3, n);
    fball_add(res, v + 0, v + 2, n);

    for (i = 0; i < 3; i++)
        fball_clear(v + i);
}

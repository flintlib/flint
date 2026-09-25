/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mp_real.h"
#include "impl.h"

/* log 10 = log 2 + log 5, the two from the three-prime Machin-type set
   (46 atanh(1/31) + 34 atanh(1/49) + 20 atanh(1/161) in the plain
   arctangent form, each term from Zuniga's series where that is
   faster), which gives log 2, log 3 and log 5 from the same three
   series and runs them on the available threads. */
void
_mp_real_const_log10_compute(mp_real_t res, slong n)
{
    mp_real_struct v[3];
    slong i;

    for (i = 0; i < 3; i++)
        mp_real_init(v + i);

    _mp_real_log_primes_vec(v, 3, n);
    mp_real_add(res, v + 0, v + 2, n);

    for (i = 0; i < 3; i++)
        mp_real_clear(v + i);
}

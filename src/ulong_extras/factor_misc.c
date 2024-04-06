/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong n_factor_evaluate(const n_factor_t * fac)
{
    ulong ret = 1;
    slong ix;

    for (ix = 0; ix < fac->num; ix++)
    {
        ulong t0, t1;
        t0 = _n_pow_check(fac->p[ix], fac->exp[ix]);
        umul_ppmm(t1, ret, ret, t0);
        if (t1 != 0)
            return 0;
    }

    return ret;
}

/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

slong
acb_dirichlet_platt_hardy_z_zeros(
        arb_ptr res, const fmpz_t n, slong len, slong prec)
{
    if (len <= 0 || fmpz_sizeinbase(n, 10) < 5)
    {
        return 0;
    }
    else if (fmpz_sgn(n) < 1)
    {
        flint_throw(FLINT_ERROR, "Nonpositive indices of Hardy Z zeros are not supported.\n");
    }
    else
    {
        slong r, s;
        fmpz_t k;
        fmpz_init(k);
        fmpz_set(k, n);
        for (s = 0; s < len; s += r)
        {
            r = acb_dirichlet_platt_local_hardy_z_zeros(res + s, k, len - s, prec);
            if (!r)
                break;
            fmpz_add_si(k, k, r);
        }
        return s;
    }
    return 0;
}

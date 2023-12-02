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
acb_dirichlet_platt_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)
{
    if (len <= 0 || fmpz_sizeinbase(n, 10) < 5)
    {
        return 0;
    }
    else if (fmpz_sgn(n) < 1)
    {
        flint_throw(FLINT_ERROR, "nonpositive indices of zeta zeros are not supported\n");
    }
    else
    {
        slong i, found;
        arb_ptr p;
        p = _arb_vec_init(len);
        found = acb_dirichlet_platt_hardy_z_zeros(p, n, len, prec);
        for (i = 0; i < found; i++)
        {
            acb_set_d(res + i, 0.5);
            arb_set(acb_imagref(res + i), p + i);
        }
        _arb_vec_clear(p, len);
        return found;
    }
    return 0;
}

/*
    Copyright (C) 2010 Juan Arias de Reyna
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)
{
    if (len <= 0)
    {
        return;
    }
    else if (fmpz_sgn(n) < 1)
    {
        flint_throw(FLINT_ERROR, "nonpositive indices of zeta zeros are not supported\n");
    }
    else
    {
        slong i;
        arb_ptr p;
        p = _arb_vec_init(len);
        acb_dirichlet_hardy_z_zeros(p, n, len, prec);
        for (i = 0; i < len; i++)
        {
            acb_set_d(res + i, 0.5);
            arb_set(acb_imagref(res + i), p + i);
        }
        _arb_vec_clear(p, len);
    }
}

/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_agm_mul_all(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr v;
    ulong a, b;

    v = _acb_vec_init(n);

    for (a = 0; a < n; a++)
    {
        _acb_vec_set(v, a2, n);
        for (b = 0; b < n; b++)
        {
            if (acb_theta_char_dot(a, b, g) % 2 == 1)
            {
                acb_neg(&v[b], &v[b]);
            }
        }
        acb_theta_agm_mul(v, a1, v, g, prec);
        for (b = 0; b < n; b++)
        {
            acb_set(&res[b * n + a], &v[b]);
        }
    }
    _acb_vec_scalar_mul_2exp_si(res, res, n * n, g);

    _acb_vec_clear(v, n);
}

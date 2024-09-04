/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static void
acb_theta_agm_hadamard(acb_ptr res, acb_srcptr a, slong g, slong prec)
{
    acb_ptr v;
    slong half;

    if (g == 0)
    {
        acb_set(&res[0], &a[0]);
    }
    else
    {
        half = 1 << (g - 1);
        v = _acb_vec_init(1 << g);

        acb_theta_agm_hadamard(v, a, g - 1, prec);
        acb_theta_agm_hadamard(v + half, a + half, g - 1, prec);
        _acb_vec_add(res, v, v + half, half, prec);
        _acb_vec_sub(res + half, v, v + half, half, prec);

        _acb_vec_clear(v, 1 << g);
    }
}

void
acb_theta_agm_mul(acb_ptr res, acb_srcptr a1, acb_srcptr a2, slong g, slong prec)
{
    acb_ptr v;
    slong n = 1 << g;
    slong k;

    v = _acb_vec_init(2 * n);

    acb_theta_agm_hadamard(v, a1, g, prec);

    if (a1 == a2)
    {
        for (k = 0; k < n; k++)
        {
            acb_sqr(&v[k], &v[k], prec);
        }
    }
    else
    {
        acb_theta_agm_hadamard(v + n, a2, g, prec);
        for (k = 0; k < n; k++)
        {
            acb_mul(&v[k], &v[k], &v[k + n], prec);
        }
    }

    acb_theta_agm_hadamard(res, v, g, prec);
    _acb_vec_scalar_mul_2exp_si(res, res, n, -2 * g);

    _acb_vec_clear(v, 2 * n);
}

/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_agm_mul(acb_ptr r, acb_srcptr a1, acb_srcptr a2, slong g, slong prec)
{
    acb_ptr v;
    slong n = 1 << g;
    slong k;

    v = _acb_vec_init(2 * n);

    acb_theta_agm_hadamard(v, a1, g, prec);
    acb_theta_agm_hadamard(v + n, a2, g, prec);
    for (k = 0; k < n; k++)
    {
        acb_mul(&v[k], &v[k], &v[k + n], prec);
    }
    acb_theta_agm_hadamard(r, v, g, prec);
    _acb_vec_scalar_mul_2exp_si(r, r, n, -2 * g);

    _acb_vec_clear(v, 2 * n);
}
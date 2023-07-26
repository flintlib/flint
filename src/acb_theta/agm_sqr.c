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
acb_theta_agm_sqr(acb_ptr r, acb_srcptr a, slong g, slong prec)
{
    acb_ptr v;
    slong n = 1 << g;

    v = _acb_vec_init(n);

    acb_theta_agm_hadamard(v, a, g, prec);
    _acb_vec_sqr(v, v, n, prec);
    acb_theta_agm_hadamard(r, v, g, prec);
    _acb_vec_scalar_mul_2exp_si(r, r, n, -2 * g);

    _acb_vec_clear(v, n);
}

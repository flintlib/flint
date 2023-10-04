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
acb_theta_ql_step_2(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;

    aux = _acb_vec_init(3 * n);

    /* Duplication using square roots for z + t and z + 2t */
    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    _acb_vec_set(res, aux, 3 * n);

    _acb_vec_clear(aux, 3 * n);
}

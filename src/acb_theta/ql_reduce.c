/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong acb_theta_ql_reduce(acb_ptr x, acb_t c, arb_t u, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t cho;
    arb_ptr offset;
    arf_t R2, eps;
    arb_t bound, t;
    slong d;

    arb_mat_init(cho, g, g);
    offset = _arb_vec_init(g);
    arf_init(R2);
    arf_init(eps);
    arb_init(bound);
    arb_init(t);

    acb_theta_eld_cho(cho, tau, prec);
    acb_theta_naive_radius(R2, eps, cho, 0, prec);
    acb_theta_naive_reduce(offset, x, c, u, z, 1, tau, cho, prec);
    arb_mul_arf(u, u, eps, prec);

    arb_set_arf(bound, R2);
    arb_mul_2exp_si(bound, bound, 2);

    for (d = g; d > 0; d--)
    {
        arb_sqr(t, arb_mat_entry(cho, d - 1, d - 1), prec);
        if (!arb_gt(t, bound))
        {
            break;
        }
    }

    arb_mat_clear(cho);
    _arb_vec_clear(offset, g);
    arf_clear(R2);
    arf_clear(eps);
    arb_clear(bound);
    arb_clear(t);
    return d;
}

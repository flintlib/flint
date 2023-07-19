/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_ql_dist_rec(arb_t x, arb_srcptr offset, const arb_mat_t cho, slong d, slong prec);

static void
acb_theta_ql_dist_fixed_coord(arb_t x, arb_srcptr offset, slong n,
    const arb_mat_t cho, slong d, slong prec)
{    
    arb_ptr new_offset;
    arb_t c;
    slong k;

    new_offset = _arb_vec_init(d - 1);
    arb_init(c);

    arb_set(c, &offset[d - 1]);
    arb_div(c, c, arb_mat_entry(cho, d - 1, d - 1), prec);
    arb_sub_si(c, c, n, prec);

    for (k = 0; k < d - 1; k++)
    {
        arb_mul(&new_offset[k], arb_mat_entry(cho, k, d - 1), c, prec);
    }
    _arb_vec_sub(new_offset, offset, new_offset, d - 1, prec);
    
    acb_theta_ql_dist_rec(x, new_offset, cho, d - 1, prec);
    arb_mul(c, c, arb_mat_entry(cho, d - 1, d - 1), prec);
    arb_sqr(c, c, prec);
    arb_add(x, x, c, prec);
    
    _arb_vec_clear(new_offset, d - 1);
    arb_clear(c);
}

static void
acb_theta_ql_dist_rec(arb_t x, arb_srcptr offset, const arb_mat_t cho, slong d, slong prec)
{
    arb_t c, y;
    arf_t rad;
    slong min, mid, max;
    slong k;

    if (d == 0)
    {
        arb_zero(x);
        return;
    }

    arb_init(c);
    arb_init(y);
    arf_init(rad);
    
    arb_set(c, &offset[d - 1]);
    arb_div(c, c, arb_mat_entry(cho, d - 1, d - 1), prec);
    mid = arf_get_si(arb_midref(c), ARF_RND_NEAR);

    acb_theta_ql_dist_fixed_coord(y, offset, mid, cho, d, prec);
    arb_get_ubound_arf(rad, y, prec);

    /* eld_interval uses even integers only */
    arb_mul_2exp_si(c, c, 1);
    arf_mul_2exp_si(rad, rad, 1);
    acb_theta_eld_interval(&min, &mid, &max, c, rad, 0, prec);

    arb_set(x, y);
    for (k = min/2; k <= max/2; k++)
    {
        acb_theta_ql_dist_fixed_coord(y, offset, k, cho, d, prec);
        arb_min(x, x, y, prec);
    }

    arb_clear(c);
    arb_clear(y);
    arf_clear(rad);
}

void
acb_theta_ql_dist(arb_t x, arb_srcptr offset, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    acb_theta_ql_dist_rec(x, offset, cho, g, prec);
}

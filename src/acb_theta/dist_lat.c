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
acb_theta_dist_ubound(arf_t u, arb_srcptr v, const arb_mat_t cho, slong prec)
{
    slong g = acb_mat_nrows(cho);
    slong nb = 1 << g;
    arb_mat_t m;
    arb_ptr x;
    slong* approx;
    slong* pt;
    arb_t d2;
    arf_t b;
    slong j, k;

    arb_mat_init(m, g, g);
    x = _arb_vec_init(g);
    approx = flint_malloc(2 * g * sizeof(slong));
    pt = flint_malloc(g * sizeof(slong));
    arb_init(d2);
    arf_init(b);

    arb_mat_inv(m, cho, prec);
    arb_mat_vector_mul_col(x, m, v, prec);
    for (k = 0; k < g; k++)
    {
        approx[2 * k] = - arf_get_si(arb_midref(&x[k]), ARF_RND_FLOOR);
        approx[2 * k + 1] = - arf_get_si(arb_midref(&x[k]), ARF_RND_CEIL);
    }

    arf_pos_inf(u);
    for (k = 0; k < nb; k++)
    {
        for (j = 0; j < g; j++)
        {
            if (k & (1 << j))
            {
                pt[j] = approx[2 * j];
            }
            else
            {
                pt[j] = approx[2 * j + 1];
            }
        }
        acb_theta_dist_pt(d2, v, cho, pt, prec);
        arb_get_ubound_arf(b, d2, prec);
        arf_min(u, u, b);
    }

    arb_mat_clear(m);
    _arb_vec_clear(x, g);
    flint_free(approx);
    flint_free(pt);
    arb_clear(d2);
    arf_clear(b);
}

void
acb_theta_dist_lat(arb_t d2, arb_srcptr v, const arb_mat_t cho, slong prec)
{
    slong g = arb_mat_nrows(cho);
    acb_theta_eld_t E;
    slong nb;
    slong* pts;
    arf_t u;
    arb_t x;
    slong k;

    acb_theta_eld_init(E, g, g);
    arf_init(u);
    arb_init(x);

    acb_theta_dist_ubound(u, v, cho, prec);
    acb_theta_eld_fill(E, cho, u, v, prec);
    nb = acb_theta_eld_nb_pts(E);

    pts = flint_malloc(nb * g * sizeof(slong));
    acb_theta_eld_points(pts, E);

    arb_pos_inf(d2);
    for (k = 0; k < nb; k++)
    {
        acb_theta_dist_pt(x, v, cho, pts + k * g, prec);
        arb_min(d2, d2, x, prec);
    }

    acb_theta_eld_clear(E);
    arf_clear(u);
    arb_clear(x);
    flint_free(pts);
}

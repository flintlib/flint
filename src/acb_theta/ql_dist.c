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
acb_theta_ql_dist(arb_t d2, arb_srcptr v, const arb_mat_t cho, slong prec)
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

    acb_theta_ql_dist_ubound(u, v, cho, prec);
    acb_theta_eld_fill(E, cho, u, v, prec);
    nb = acb_theta_eld_nb_pts(E);

    pts = flint_malloc(nb * g * sizeof(slong));
    acb_theta_eld_points(pts, E);

    arb_pos_inf(d2);
    for (k = 0; k < nb; k++)
    {
        acb_theta_ql_dist_pt(x, v, cho, pts + k * g, prec);
        arb_min(d2, d2, x, prec);
    }

    acb_theta_eld_clear(E);
    arf_clear(u);
    arb_clear(x);
    flint_free(pts);    
}

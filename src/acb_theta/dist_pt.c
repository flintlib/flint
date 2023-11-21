/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_theta.h"

void
acb_theta_dist_pt(arb_t d, arb_srcptr v, const arb_mat_t C, const slong * n, slong prec)
{
    slong g = arb_mat_nrows(C);
    arb_ptr w;
    slong k;

    w = _arb_vec_init(g);

    for (k = 0; k < g; k++)
    {
        arb_set_si(&w[k], n[k]);
    }
    arb_mat_vector_mul_col(w, C, w, prec);
    _arb_vec_add(w, w, v, g, prec);
    arb_dot(d, NULL, 0, w, 1, w, 1, g, prec);

    _arb_vec_clear(w, g);
}

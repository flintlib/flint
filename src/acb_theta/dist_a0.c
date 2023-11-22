/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_dist_a0(arb_ptr d, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t Yinv, C;
    arb_ptr v, w;
    ulong a;

    arb_mat_init(Yinv, g, g);
    arb_mat_init(C, g, g);
    v = _arb_vec_init(g);
    w = _arb_vec_init(g);

    acb_siegel_yinv(Yinv, tau, prec);
    acb_siegel_cho(C, tau, prec);

    _acb_vec_get_imag(v, z, g);
    arb_mat_vector_mul_col(v, Yinv, v, prec);

    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_arb(w, a, g);
        _arb_vec_add(w, v, w, g, prec);
        arb_mat_vector_mul_col(w, C, w, prec);
        acb_theta_dist_lat(&d[a], w, C, prec);
    }

    arb_mat_clear(Yinv);
    arb_mat_clear(C);
    _arb_vec_clear(v, g);
    _arb_vec_clear(w, g);
}

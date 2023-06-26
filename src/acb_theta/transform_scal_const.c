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
acb_theta_transform_scal_const(acb_t scal, const acb_mat_t tau,
                               const fmpz_mat_t mat, slong k2, slong prec)
{
    slong g = acb_mat_nrows(tau);
    acb_t mu;
    acb_t det;
    acb_mat_t w;

    acb_init(mu);
    acb_init(det);
    acb_mat_init(w, g, g);

    acb_onei(mu);
    acb_pow_si(mu, mu, k2, prec);
    acb_siegel_cocycle(w, mat, tau, prec);
    acb_mat_det(det, w, prec);
    acb_mul(scal, det, mu, prec);

    acb_clear(mu);
    acb_clear(det);
    acb_mat_clear(w);
}

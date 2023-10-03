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
acb_siegel_cocycle_det(acb_t res, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = sp2gz_dim(mat);
    acb_mat_t w;

    acb_mat_init(w, g, g);

    acb_siegel_cocycle(w, mat, tau, prec);
    acb_mat_det(res, w, prec);

    acb_mat_clear(w);
}

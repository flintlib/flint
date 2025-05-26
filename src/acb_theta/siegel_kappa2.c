/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

slong
acb_siegel_kappa2(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong prec = 2;
    acb_mat_t tau;
    acb_t s;
    slong res;

    acb_mat_init(tau, g, g);
    acb_init(s);

    res = acb_siegel_kappa(s, mat, tau, 1, prec);

    acb_mat_clear(tau);
    acb_clear(s);
    return res;
}

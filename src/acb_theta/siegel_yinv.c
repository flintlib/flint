/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_yinv(arb_mat_t Yinv, const acb_mat_t tau, slong prec)
{
    int res;

    acb_mat_get_imag(Yinv, tau);
    res = arb_mat_inv(Yinv, Yinv, prec);
    if (!res)
    {
        arb_mat_indeterminate(Yinv);
    }
}

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
acb_siegel_cho(arb_mat_t C, const acb_mat_t tau, slong prec)
{
    arb_t pi;
    int res;

    arb_init(pi);
    arb_const_pi(pi, prec);

    acb_mat_get_imag(C, tau);
    arb_mat_scalar_mul_arb(C, C, pi, prec);
    res = arb_mat_cho(C, C, prec);
    arb_mat_transpose(C, C);

    if (!res)
    {
        arb_mat_indeterminate(C);
    }

    arb_clear(pi);
}

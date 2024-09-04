/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_cho_yinv(arb_mat_t cho, arb_mat_t yinv, const acb_mat_t tau, slong prec)
{
    arb_t sqrtpi;
    int res;

    arb_init(sqrtpi);
    arb_const_sqrt_pi(sqrtpi, prec);

    acb_mat_get_imag(cho, tau);
    res = arb_mat_cho(cho, cho, prec);
    if (!res)
    {
        arb_mat_indeterminate(cho);
    }
    arb_mat_inv_cho_precomp(yinv, cho, prec);

    arb_mat_transpose(cho, cho);
    arb_mat_scalar_mul_arb(cho, cho, sqrtpi, prec);

    arb_clear(sqrtpi);
}

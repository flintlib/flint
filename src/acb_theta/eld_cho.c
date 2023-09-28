/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_eld_cho(arb_mat_t cho, const acb_mat_t tau, slong prec)
{
    arb_t pi;
    int res;

    arb_init(pi);
    arb_const_pi(pi, prec);

    acb_mat_get_imag(cho, tau);
    arb_mat_scalar_mul_arb(cho, cho, pi, prec);
    res = arb_mat_cho(cho, cho, prec);
    arb_mat_transpose(cho, cho);

    if (!res)
    {
        arb_mat_indeterminate(cho);
    }

    arb_clear(pi);
}

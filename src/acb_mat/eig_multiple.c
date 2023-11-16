/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int
acb_mat_eig_multiple(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
{
    slong n;
    acb_ptr F;
    int success;

    n = arb_mat_nrows(A);
    F = _acb_vec_init(n);

    success = acb_mat_eig_simple_vdhoeven_mourrain(F, NULL, NULL, A, E_approx, R_approx, prec);

    if (!success)
        success = acb_mat_eig_multiple_rump(F, A, E_approx, R_approx, prec);

    _acb_vec_set(E, F, n);
    _acb_vec_clear(F, n);

    return success;
}

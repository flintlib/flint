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
acb_mat_eig_simple(acb_ptr E, acb_mat_t L, acb_mat_t R,
    const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
{
    return acb_mat_eig_simple_vdhoeven_mourrain(E, L, R, A, E_approx, R_approx, prec);
}

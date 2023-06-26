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
acb_siegel_reduce_imag(fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t im;
    fmpz_mat_t U;

    arb_mat_init(im, g, g);
    fmpz_mat_init(U, g, g);

    acb_mat_get_imag(im, tau);
    arb_mat_reduce(U, im, prec);

    fmpz_mat_diag_sp(mat, U);

    arb_mat_clear(im);
    fmpz_mat_clear(U);
}

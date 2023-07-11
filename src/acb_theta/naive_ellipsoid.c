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
acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr c, acb_ptr new_z,
    ulong ab, int all, slong ord, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, const arf_t eps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong eld_prec = ACB_THETA_ELD_DEFAULT_PREC;
    arb_t pi;
    arf_t R2;
    slong scl = -1;
    arb_mat_t cho;
    arb_ptr offset;
    int res;

    arb_init(pi);
    arf_init(R2);
    arb_mat_init(cho, g, g);
    offset = _arb_vec_init(g);
    
    if (all)
    {
        ab = 0;
        scl = -2;
    }

    acb_mat_get_imag(cho, tau);
    arb_const_pi(pi, prec);
    arb_mat_scalar_mul_arb(cho, cho, pi, prec);

    /* Get Cholesky for pi Y, possibly at high precision */
    res = arb_mat_cho(cho, cho, eld_prec);
    if (!res)
    {
        eld_prec = prec;
        res = arb_mat_cho(cho, cho, eld_prec);
    }
    if (!res)
    {
        flint_printf("acb_theta_naive_ellipsoid: Error ");
        flint_printf("(imaginary part is not positive definite)\n");
        fflush(stdout);
        flint_abort();
    }
    arb_mat_transpose(cho, cho);

    /* Get radius for error of at most eps */
    acb_theta_naive_radius(R2, cho, ord, eps, eld_prec);

    /* Reduce all z, set offset */
    acb_theta_naive_reduce(offset, new_z, c, z, nb_z, tau, cho, prec);

    /* Fill ellipsoid */
    arb_mat_scalar_mul_2exp_si(cho, cho, scl);
    acb_theta_eld_fill(E, cho, R2, offset, NULL, ab >> g, eld_prec);

    arb_clear(pi);
    arf_clear(R2);
    arb_mat_clear(cho);
    _arb_vec_clear(offset, g);
}

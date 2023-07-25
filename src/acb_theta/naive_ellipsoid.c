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
acb_theta_naive_ellipsoid(acb_theta_eld_t E, acb_ptr new_z, acb_ptr c, arb_ptr u,
    slong ord, acb_srcptr z, slong nb_z, const acb_mat_t tau, const arf_t eps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong eld_prec = ACB_THETA_ELD_DEFAULT_PREC;
    arb_t pi;
    arf_t R2;
    arb_mat_t cho;
    arb_ptr offset;
    int res;
    slong k;

    arb_init(pi);
    arf_init(R2);
    arb_mat_init(cho, g, g);
    offset = _arb_vec_init(g);
    
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
        acb_mat_printd(tau, 5);
        fflush(stdout);
        flint_abort();
    }
    arb_mat_transpose(cho, cho);

    /* Reduce all z, set offset and upper bounds */
    acb_theta_naive_reduce(offset, new_z, c, u, z, nb_z, tau, cho, prec);
    for (k = 0; k < nb_z; k++)
    {
        arb_mul_arf(&u[k], &u[k], eps, prec);
    }
    
    /* Get radius for error of at most eps and fill ellipsoid */
    acb_theta_naive_radius(R2, cho, ord, eps, eld_prec);
    acb_theta_eld_fill(E, cho, R2, offset, eld_prec);

    arb_clear(pi);
    arf_clear(R2);
    arb_mat_clear(cho);
    _arb_vec_clear(offset, g);
}
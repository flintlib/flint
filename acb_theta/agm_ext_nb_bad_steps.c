/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_agm_ext_nb_bad_steps(acb_srcptr z, const acb_mat_t tau, slong prec)
{
    arb_mat_t im;
    arb_t lambda;
    arb_t lambda0;
    arb_t norm;
    arb_t temp;
    arf_t up;
    fmpz_t e;
    slong g = acb_mat_nrows(tau);
    slong k;
    slong res;

    arb_mat_init(im, g, g);
    arb_init(lambda);
    arb_init(lambda0);
    arb_init(norm);
    arb_init(temp);
    arf_init(up);
    fmpz_init(e);

    /* Get lambda = smallest eigenvalue of Pi Im(tau) */
    acb_mat_get_imag(im, tau);
    arb_mat_pos_lambda(lambda, im, prec);
    arb_const_pi(lambda0, prec);
    arb_mul(lambda, lambda, lambda0, prec);

    /* Set lambda0 such that 3g exp(-lambda0) = 1/50 */
    arb_one(lambda0);
    arb_div_si(lambda0, lambda0, 150 * g, prec);
    arb_log(lambda0, lambda0, prec);
    arb_neg(lambda0, lambda0);

    /* Max lambda0 with 4*norm(Im(z)) */
    arb_zero(norm);
    for (k = 0; k < g; k++)
    {
        arb_set(temp, acb_imagref(&z[k]));
        arb_sqr(temp, temp, prec);
        arb_add(norm, norm, temp, prec);
    }
    arb_sqrt(norm, norm, prec);
    arb_mul_2exp_si(norm, norm, 2);
    arb_max(lambda0, lambda0, norm, prec);

    /* Compute n, minimal s.t. 2^n lambda > 2lambda0 */
    arb_div(lambda, lambda0, lambda, prec);
    arb_get_ubound_arf(up, lambda, prec);

    if (!arf_is_finite(up))
    {
        flint_printf("agm_ext_nb_bad_steps: Error (infinite value)\n");
        fflush(stdout);
        flint_abort();
    }

    arf_frexp(up, e, up);
    res = fmpz_get_si(e) + 1;

    arb_mat_clear(im);
    arb_clear(lambda);
    arb_clear(lambda0);
    arb_clear(norm);
    arb_clear(temp);
    arf_clear(up);
    fmpz_clear(e);
    return res;
}

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

/* Compute c0, c1, c2 such that |theta_{a,b}(z,tau)| on a ball of radius rho
   around z is bounded above by c0 exp((c1 + c2 rho)^2) */

static void
acb_theta_jet_ql_ci(arb_t c0, arb_t c1, arb_t c2, acb_srcptr z, const acb_mat_t tau)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong g = acb_mat_nrows(tau);
    arb_mat_t Yinv;
    arb_mat_t cho;
    arb_ptr y, w;
    arb_t t, s;
    slong k, j;

    arb_mat_init(Yinv, g, g);
    arb_mat_init(cho, g, g);
    y = _arb_vec_init(g);
    w = _arb_vec_init(g);
    arb_init(t);
    arb_init(s);

    _acb_vec_get_imag(y, z, g);
    acb_siegel_yinv(Yinv, tau, lp);
    acb_siegel_cho(cho, tau, lp);

    /* c0 is 2^g \prod_{i=1}^g (1 + 2/\sqrt{\gamma_i}) */
    arb_one(c0);
    arb_mul_2exp_si(c0, c0, g);
    for (k = 0; k < g; k++)
    {
        arb_mul_2exp_si(t, arb_mat_entry(cho, k, k), 1);
        arb_add_si(t, t, 1, lp);
        arb_mul(c0, c0, t, lp);
    }

    /* c1 is sqrt(\pi y Y^{-1} y) */
    arb_const_pi(t, lp);
    arb_mat_scalar_mul_arb(Yinv, Yinv, t, lp);
    arb_mat_vector_mul_col(w, Yinv, y, lp);
    arb_dot(c1, NULL, 0, y, 1, w, 1, g, lp);
    arb_nonnegative_part(c1, c1);
    arb_sqrt(c1, c1, lp);

    /* c2 is sqrt(max of \pi x Y^{-1} x where |x| \leq 1) */
    arb_zero(c2);
    arb_mat_cho(cho, Yinv, lp);
    arb_mat_transpose(cho, cho);
    for (k = 0; k < g; k++)
    {
        arb_zero(s);
        for (j = k; j < g; j++)
        {
            arb_abs(t, arb_mat_entry(cho, k, j));
            arb_add(s, s, t, lp);
        }
        arb_sqr(s, s, lp);
        arb_add(c2, c2, s, lp);
    }
    arb_nonnegative_part(c2, c2);
    arb_sqrt(c2, c2, lp);

    arb_mat_clear(Yinv);
    arb_mat_clear(cho);
    _arb_vec_clear(y, g);
    _arb_vec_clear(w, g);
    arb_clear(t);
    arb_clear(s);
}

/* Pick rho and c such that |theta_{a,b}(z,tau)| on a ball of radius rho around
   z is bounded above by c, and is a good choice to bound derivatives up to
   order ord */

void
acb_theta_jet_ql_bounds(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong ord)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong b = ord + 1;
    arb_t t, c0, c1, c2;
    arf_t x;

    arb_init(t);
    arb_init(c0);
    arb_init(c1);
    arb_init(c2);
    arf_init(x);

    /* Get ci */
    acb_theta_jet_ql_ci(c0, c1, c2, z, tau);

    /* Set rho to positive root of 2 c_2 rho (c_1 + c_2 rho) = 2 b - 1 */
    arb_mul(t, c1, c2, lp);
    arb_mul_2exp_si(t, t, 1);
    arb_sqr(rho, t, lp);
    arb_sqr(t, c2, lp);
    arb_mul_si(t, t, 8 * (2 * b - 1), lp);
    arb_add(rho, rho, t, lp);
    arb_sqrt(rho, rho, lp);
    arb_mul(t, c1, c2, lp);
    arb_submul_si(rho, t, 2, lp);
    arb_sqr(t, c2, lp);
    arb_mul_2exp_si(t, t, 2);
    arb_div(rho, rho, t, lp);

    /* Set c to corresponding bound */
    arb_mul(c, c2, rho, lp);
    arb_add(c, c, c1, lp);
    arb_sqr(c, c, lp);
    arb_exp(c, c, lp);
    arb_mul(c, c, c0, lp);

    arb_clear(t);
    arb_clear(c0);
    arb_clear(c1);
    arb_clear(c2);
    arf_clear(x);
}

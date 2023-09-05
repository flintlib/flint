/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Compute c0, c1, c2 such that |theta_{a,b}(z,tau)| on a ball of radius rho
   around z is bounded above by c0 exp((c1 + c2 rho)^2) */

static void
acb_theta_jet_bounds_ci(arb_t c0, arb_t c1, arb_t c2, acb_srcptr z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t Yinv;
    arb_mat_t cho;
    arb_ptr y;
    arb_t t, s;
    slong k, j;

    arb_mat_init(Yinv, g, g);
    arb_mat_init(cho, g, g);
    y = _arb_vec_init(g);
    arb_init(t);
    arb_init(s);

    acb_mat_get_imag(Yinv, tau);
    _acb_vec_get_imag(y, z, g);
    arb_mat_inv(Yinv, Yinv, prec);
    acb_theta_eld_cho(cho, tau, prec);

    /* c0 is 2^g \prod_{i=1}^g (1 + 2/\sqrt{\gamma_i}) */
    arb_one(c0);
    arb_mul_2exp_si(c0, c0, g);
    for (k = 0; k < g; k++)
    {
        arb_mul_2exp_si(t, arb_mat_entry(cho, k, k), 1);
        arb_add_si(t, t, 1, prec);
        arb_mul(c0, c0, t, prec);
    }

    /* c1 is sqrt(\pi y Y^{-1} y) */
    arb_const_pi(t, prec);
    arb_mat_scalar_mul_arb(Yinv, Yinv, t, prec);
    arb_mat_bilinear_form(c1, Yinv, y, y, prec);
    arb_sqrt(c1, c1, prec);

    /* c2 is sqrt(max of \pi x Y^{-1} x where |x| \leq 1) */
    arb_zero(c2);
    arb_mat_cho(cho, Yinv, prec);
    arb_mat_transpose(cho, cho);
    for (k = 0; k < g; k++)
    {
        arb_zero(s);
        for (j = k; j < g; j++)
        {
            arb_abs(t, arb_mat_entry(cho, k, j));
            arb_add(s, s, t, prec);
        }
        arb_sqr(s, s, prec);
        arb_add(c2, c2, s, prec);
    }
    arb_sqrt(c2, c2, prec);

    arb_mat_clear(Yinv);
    arb_mat_clear(cho);
    _arb_vec_clear(y, g);
    arb_clear(t);
    arb_clear(s);
}

/* Pick rho and c such that |theta_{a,b}(z,tau)| on a ball of radius rho around
   z is bounded above by c, and is a good choice to bound derivatives up to
   order ord */

void
acb_theta_jet_bounds(arb_t eps, arb_t c, arb_t rho, acb_srcptr z,
    const acb_mat_t tau, slong ord, slong hprec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t t, c0, c1, c2;
    arf_t x;

    arb_init(t);
    arb_init(c0);
    arb_init(c1);
    arb_init(c2);
    arf_init(x);

    /* Get ci */
    acb_theta_jet_bounds_ci(c0, c1, c2, z, tau, prec);

    /* Set rho to positive root of 2 c_2 rho (c_1 + c_2 rho) = 2 ord */
    arb_mul(t, c1, c2, prec);
    arb_mul_2exp_si(t, t, 1);
    arb_sqr(rho, t, prec);
    arb_sqr(t, c2, prec);
    arb_mul_si(t, t, 64 * ord, prec);
    arb_add(rho, rho, t, prec);
    arb_sqrt(rho, rho, prec);
    arb_mul(t, c1, c2, prec);
    arb_submul_si(rho, t, 2, prec);
    arb_sqr(t, c2, prec);
    arb_mul_2exp_si(t, t, 2);
    arb_div(rho, rho, t, prec);

    /* Set c to corresponding bound */
    arb_mul(c, c2, rho, prec);
    arb_add(c, c, c1, prec);
    arb_sqr(c, c, prec);
    arb_exp(c, c, prec);
    arb_mul(c, c, c0, prec);

    /* Set eps to minimum of rho/g^(1/ord) and (2^(-hprec)/cg)^{1/n} rho^2 */
    arb_set_si(eps, g);
    arb_root_ui(eps, eps, ord, prec);
    arb_mul_si(t, c, g, prec);
    arb_inv(t, t, prec);
    arb_mul_2exp_si(t, t, -hprec);
    arb_root_ui(t, t, ord, prec);
    arb_mul(t, t, rho, prec);
    arb_min(eps, eps, t, prec);
    arb_mul(eps, eps, rho, prec);
    arb_get_lbound_arf(x, eps, prec);
    arb_set_arf(eps, x); /* eps is exact */

    /* Set c to c * 2^(- hprec + 1) */
    arb_mul_2exp_si(c, c, -hprec + 1);

    arb_clear(t);
    arb_clear(c0);
    arb_clear(c1);
    arb_clear(c2);
    arf_clear(x);
}

/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_hypgeom.h"

/* c1 = 4*zeta(sigma)/(2*sigma-1) * pi^((-3-2*sigma)/4) * exp(...) */
static void
_pre_i_c1(arb_t res, slong A, const arb_t H, slong sigma, slong prec)
{
    arb_t pi, x1, x2, x3, x4;

    arb_init(pi);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);

    arb_const_pi(pi, prec);

    /* x1 = 4*zeta(sigma) / (2*sigma - 1) */
    arb_set_si(x1, sigma);
    arb_zeta(x1, x1, prec);
    arb_mul_2exp_si(x1, x1, 2);
    arb_div_si(x1, x1, 2*sigma - 1, prec);

    /* x2 = pi^((-3-2*sigma)/4) */
    arb_inv(x2, pi, prec);
    arb_root_ui(x2, x2, 4, prec);
    arb_pow_ui(x2, x2, 3 + 2*sigma, prec);

    /* x3 = (2*sigma-1)^2 / (8*H^2) */
    arb_set_si(x3, 2*sigma - 1);
    arb_div(x3, x3, H, prec);
    arb_sqr(x3, x3, prec);
    arb_mul_2exp_si(x3, x3, -3);

    /* x4 = pi*A*(2*sigma-1)/2 */
    arb_mul_si(x4, pi, A*(2*sigma-1), prec);
    arb_mul_2exp_si(x4, x4, -1);

    /* res = x1*x2*exp(x3 - x4) */
    arb_sub(res, x3, x4, prec);
    arb_exp(res, res, prec);
    arb_mul(res, res, x1, prec);
    arb_mul(res, res, x2, prec);

    arb_clear(pi);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
}

/* c2 = 8*pi^(1/4) */
static void
_pre_i_c2(arb_t res, slong prec)
{
    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);
    arb_root_ui(res, pi, 4, prec);
    arb_mul_2exp_si(res, res, 3);
    arb_clear(pi);
}

void
acb_dirichlet_platt_i_precomp_init(acb_dirichlet_platt_i_precomp_t pre,
        slong A, const arb_t H, slong sigma, slong prec)
{
    arb_init(&pre->c1);
    arb_init(&pre->c2);
    _pre_i_c1(&pre->c1, A, H, sigma, prec);
    _pre_i_c2(&pre->c2, prec);
}

void
acb_dirichlet_platt_i_precomp_clear(acb_dirichlet_platt_i_precomp_t pre)
{
    arb_clear(&pre->c1);
    arb_clear(&pre->c2);
}

void
acb_dirichlet_platt_i_bound_precomp(arb_t res,
        const acb_dirichlet_platt_i_precomp_t pre_i,
        const acb_dirichlet_platt_c_precomp_t pre_c,
        const arb_t t0, slong A, const arb_t H,
        slong sigma, slong prec)
{
    arb_t pi, cbound, x1, x2, x3, x4, x5;

    arb_init(pi);
    arb_init(cbound);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_const_pi(pi, prec);

    /* x1 = (1 - 4*t0^2)/(8*H^2) */
    arb_sqr(x1, t0, prec);
    arb_mul_2exp_si(x1, x1, 2);
    arb_sub_ui(x1, x1, 1, prec);
    arb_neg(x1, x1);
    arb_div(x1, x1, H, prec);
    arb_div(x1, x1, H, prec);
    arb_mul_2exp_si(x1, x1, -3);

    /* x2 = pi*A/2 */
    arb_mul_si(x2, pi, A, prec);
    arb_mul_2exp_si(x2, x2, -1);

    /* x3 = exp(x1 - x2) */
    arb_sub(x3, x1, x2, prec);
    arb_exp(x3, x3, prec);

    acb_dirichlet_platt_c_bound_precomp(
            cbound, pre_c, sigma, t0, H, 0, prec);

    /* res = c1*cbound + c2*x3 */
    arb_mul(x4, &pre_i->c1, cbound, prec);
    arb_mul(x5, &pre_i->c2, x3, prec);
    arb_add(res, x4, x5, prec);

    arb_clear(pi);
    arb_clear(cbound);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

void
acb_dirichlet_platt_i_bound(arb_t res, const arb_t t0, slong A, const arb_t H,
        slong sigma, slong prec)
{
    acb_dirichlet_platt_c_precomp_t pre_c;
    acb_dirichlet_platt_i_precomp_t pre_i;

    acb_dirichlet_platt_c_precomp_init(pre_c, sigma, H, 0, prec);
    acb_dirichlet_platt_i_precomp_init(pre_i, A, H, sigma, prec);

    acb_dirichlet_platt_i_bound_precomp(
            res, pre_i, pre_c, t0, A, H, sigma, prec);

    acb_dirichlet_platt_c_precomp_clear(pre_c);
    acb_dirichlet_platt_i_precomp_clear(pre_i);
}

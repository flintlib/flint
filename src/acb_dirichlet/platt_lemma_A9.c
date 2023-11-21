/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static void
_arb_pow_d(arb_t z, const arb_t x, double d, slong prec)
{
    arb_t u;
    arb_init(u);
    arb_set_d(u, d);
    arb_pow(z, x, u, prec);
    arb_clear(u);
}

static void
_platt_lemma_A9_a(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong A, slong prec)
{
    arb_t a, pi, C;
    arb_t y1, y2, y3, y4;
    arb_t z1, z2, z3;

    arb_init(a);
    arb_init(pi);
    arb_init(C);
    arb_init(y1);
    arb_init(y2);
    arb_init(y3);
    arb_init(y4);
    arb_init(z1);
    arb_init(z2);
    arb_init(z3);

    arb_const_pi(pi, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul_si(a, a, A, prec);

    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);

    arb_set_si(y2, 2*sigma - 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);

    arb_mul_2exp_si(y3, a, -1);

    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);

    acb_dirichlet_platt_c_bound(C, sigma, t0, h, 0, prec);

    arb_zeta_ui(z1, (ulong) sigma, prec);
    arb_mul_2exp_si(z1, z1, 1);

    arb_set_si(z2, 1-2*sigma);
    arb_mul_2exp_si(z2, z2, -2);
    arb_pow(z2, pi, z2, prec);

    arb_sub(z3, y2, y3, prec);
    arb_exp(z3, z3, prec);

    arb_mul(out, z1, z2, prec);
    arb_mul(out, out, z3, prec);
    arb_mul(out, out, C, prec);
    arb_mul(out, out, y1, prec);

    arb_clear(a);
    arb_clear(pi);
    arb_clear(C);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(y3);
    arb_clear(y4);
    arb_clear(z1);
    arb_clear(z2);
    arb_clear(z3);
}


static void
_platt_lemma_A9_b(arb_t out,
        const arb_t t0, const arb_t h, slong A, slong prec)
{
    arb_t pi;
    arb_t x1, x2, x3, x4, x5;

    arb_init(pi);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_const_pi(pi, prec);

    _arb_pow_d(x1, pi, 1.25, prec);
    arb_mul_2exp_si(x1, x1, 2);

    arb_sqr(x2, t0, prec);
    arb_mul_2exp_si(x2, x2, 2);
    arb_sub_ui(x2, x2, 1, prec);
    arb_neg(x2, x2);
    arb_div(x2, x2, h, prec);
    arb_div(x2, x2, h, prec);
    arb_mul_2exp_si(x2, x2, -3);

    arb_mul_si(x3, pi, A, prec);
    arb_mul_2exp_si(x3, x3, -1);

    arb_mul_si(x4, pi, A, prec);
    arb_inv(x4, x4, prec);
    arb_add_ui(x4, x4, 1, prec);

    arb_sub(x5, x2, x3, prec);
    arb_exp(x5, x5, prec);

    arb_mul(out, x1, x4, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(pi);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}


void
acb_dirichlet_platt_lemma_A9(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong A, slong prec)
{
    arb_t a, b;

    if (sigma % 2 == 0 || sigma < 3)
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(a);
    arb_init(b);

    _platt_lemma_A9_a(a, sigma, t0, h, A, prec);
    _platt_lemma_A9_b(b, t0, h, A, prec);

    arb_add(out, a, b, prec);

    arb_clear(a);
    arb_clear(b);
}

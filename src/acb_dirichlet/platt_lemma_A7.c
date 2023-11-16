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
_platt_lemma_A7_S(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong A, slong prec)
{
    slong l;
    arb_t total, summand;
    arb_t pi, half;
    arb_t a;
    arb_t l_factorial, kd2, t02;
    arb_t x1, x2, x3, x4, x5;

    arb_init(total);
    arb_init(summand);
    arb_init(pi);
    arb_init(half);
    arb_init(a);
    arb_init(l_factorial);
    arb_init(kd2);
    arb_init(t02);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_one(half);
    arb_mul_2exp_si(half, half, -1);
    arb_const_pi(pi, prec);
    arb_one(l_factorial);
    arb_set_si(kd2, k);
    arb_mul_2exp_si(kd2, kd2, -1);
    arb_sqr(t02, t0, prec);

    for (l=0; l<=(sigma-1)/2; l++)
    {
        if (l > 1)
        {
            arb_mul_si(l_factorial, l_factorial, l, prec);
        }

        arb_mul_si(a, pi, 4*l+1, prec);
        arb_mul_si(a, a, A, prec);

        arb_inv(x1, a, prec);
        arb_add_ui(x1, x1, 1, prec);

        arb_add_si(x2, half, 2*l, prec);
        arb_sqr(x2, x2, prec);
        arb_add(x2, x2, t02, prec);
        arb_pow(x2, x2, kd2, prec);
        arb_div(x2, x2, l_factorial, prec);

        arb_set_si(x3, 4*l + 1);
        arb_div(x3, x3, h, prec);
        arb_sqr(x3, x3, prec);
        arb_mul_2exp_si(x3, x3, -3);

        arb_mul_2exp_si(x4, a, -1);

        arb_sub(x5, x3, x4, prec);
        arb_exp(x5, x5, prec);

        arb_mul(summand, x1, x2, prec);
        arb_mul(summand, summand, x5, prec);

        arb_add(total, total, summand, prec);
    }

    arb_set(out, total);

    arb_clear(total);
    arb_clear(summand);
    arb_clear(pi);
    arb_clear(half);
    arb_clear(a);
    arb_clear(l_factorial);
    arb_clear(kd2);
    arb_clear(t02);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

void
acb_dirichlet_platt_lemma_A7(arb_t out, slong sigma,
        const arb_t t0, const arb_t h, slong k, slong A, slong prec)
{
    arb_t S, C;
    arb_t pi, a;
    arb_t x1, x2;
    arb_t y1, y2, y3, y4;
    arb_t z1, z2;

    if (sigma % 2 == 0 || sigma < 3)
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(S);
    arb_init(C);
    arb_init(pi);
    arb_init(a);
    arb_init(x1);
    arb_init(x2);
    arb_init(y1);
    arb_init(y2);
    arb_init(y3);
    arb_init(y4);
    arb_init(z1);
    arb_init(z2);

    arb_const_pi(pi, prec);

    arb_pow_ui(x1, pi, (ulong) (k+1), prec);
    arb_mul_2exp_si(x1, x1, k+3);

    arb_div(x2, t0, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -1);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    _platt_lemma_A7_S(S, sigma, t0, h, k, A, prec);

    arb_mul(z1, x1, x2, prec);
    arb_mul(z1, z1, S, prec);

    arb_mul_si(a, pi, 2*sigma-1, prec);
    arb_mul_si(a, a, A, prec);

    arb_inv(y1, a, prec);
    arb_add_ui(y1, y1, 1, prec);

    arb_set_si(y2, 2*sigma + 1);
    arb_div(y2, y2, h, prec);
    arb_sqr(y2, y2, prec);
    arb_mul_2exp_si(y2, y2, -3);

    arb_mul_2exp_si(y3, a, -1);

    arb_sub(y4, y2, y3, prec);
    arb_exp(y4, y4, prec);

    acb_dirichlet_platt_c_bound(C, sigma, t0, h, k, prec);

    arb_mul(z2, y1, y4, prec);
    arb_mul(z2, z2, C, prec);
    arb_mul_2exp_si(z2, z2, 1);

    arb_add(out, z1, z2, prec);

    arb_clear(S);
    arb_clear(C);
    arb_clear(pi);
    arb_clear(a);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(y1);
    arb_clear(y2);
    arb_clear(y3);
    arb_clear(y4);
    arb_clear(z1);
    arb_clear(z2);
}

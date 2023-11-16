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

static void
_platt_lemma_A11_X(arb_t out,
        const arb_t t0, const arb_t h, slong B, const arb_t beta, slong prec)
{
    arb_t x1, x2;

    arb_init(x1);
    arb_init(x2);

    arb_set_si(x1, B);
    arb_mul_2exp_si(x1, x1, -1);
    arb_add(x1, x1, t0, prec);
    arb_pow(x1, x1, beta, prec);

    arb_set_si(x2, B);
    arb_div(x2, x2, h, prec);
    arb_sqr(x2, x2, prec);
    arb_mul_2exp_si(x2, x2, -3);
    arb_neg(x2, x2);
    arb_exp(x2, x2, prec);

    arb_mul(out, x1, x2, prec);

    arb_clear(x1);
    arb_clear(x2);
}

static void
_platt_lemma_A11_Y(arb_t out,
        const arb_t t0, const arb_t h, slong B, const arb_t beta, slong prec)
{
    arb_t x1, x2, x3, x4, x5;

    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_rsqrt_ui(x1, 2, prec);

    arb_pow(x2, t0, beta, prec);

    arb_one(x3);
    arb_mul_2exp_si(x3, x3, -1);

    arb_set_si(x4, B);
    arb_div(x4, x4, h, prec);
    arb_sqr(x4, x4, prec);
    arb_mul_2exp_si(x4, x4, -3);

    arb_hypgeom_gamma_upper(x5, x3, x4, 0, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

static void
_platt_lemma_A11_Z(arb_t out,
        const arb_t t0, const arb_t h, const arb_t beta, slong prec)
{
    arb_t two;
    arb_t x1, x2, x3, x4, x5;

    arb_init(two);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_set_ui(two, 2);

    arb_sub_ui(x1, beta, 1, prec);
    arb_mul_2exp_si(x1, x1, -1);
    arb_pow(x1, two, x1, prec);

    arb_pow(x2, h, beta, prec);

    arb_add_ui(x3, beta, 1, prec);
    arb_mul_2exp_si(x3, x3, -1);

    arb_div(x4, t0, h, prec);
    arb_sqr(x4, x4, prec);
    arb_mul_2exp_si(x4, x4, -1);

    arb_hypgeom_gamma_upper(x5, x3, x4, 0, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(two);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}

static int
_platt_lemma_A11_constraint(const arb_t t0, const arb_t h, slong B,
        const arb_t beta, slong prec)
{
    int result;

    arb_t a, b, expe;
    arb_init(a);
    arb_init(b);
    arb_init(expe);

    /* expe = exp(e) */
    arb_const_e(expe, prec);
    arb_exp(expe, expe, prec);

    /* a = beta*h^2 / t0 */
    arb_sqr(a, h, prec);
    arb_mul(a, a, beta, prec);
    arb_div(a, a, t0, prec);

    /* b = B/2 */
    arb_set_si(b, B);
    arb_mul_2exp_si(b, b, -1);

    result = arb_le(a, b) && arb_le(b, t0) && arb_gt(t0, expe);

    arb_clear(a);
    arb_clear(b);
    arb_clear(expe);

    return result;
}

void
acb_dirichlet_platt_lemma_A11(arb_t out, const arb_t t0, const arb_t h,
        slong B, slong prec)
{
    arb_t beta;
    arb_init(beta);
    acb_dirichlet_platt_beta(beta, t0, prec);

    if (_platt_lemma_A11_constraint(t0, h, B, beta, prec))
    {
        arb_t X, Y, Z;
        arb_t x1, x2;

        arb_init(X);
        arb_init(Y);
        arb_init(Z);
        arb_init(x1);
        arb_init(x2);

        _platt_lemma_A11_X(X, t0, h, B, beta, prec);
        _platt_lemma_A11_Y(Y, t0, h, B, beta, prec);
        _platt_lemma_A11_Z(Z, t0, h, beta, prec);

        arb_set_ui(x1, 2);
        arb_pow(x1, x1, beta, prec);
        arb_mul(x1, x1, h, prec);
        arb_div_si(x1, x1, B, prec);

        arb_add(x2, Y, Z, prec);
        arb_mul(x2, x2, x1, prec);
        arb_add(x2, x2, X, prec);
        arb_mul_ui(x2, x2, 6, prec);

        arb_set(out, x2);

        arb_clear(X);
        arb_clear(Y);
        arb_clear(Z);
        arb_clear(x1);
        arb_clear(x2);
    }
    else
    {
        arb_zero_pm_inf(out);
    }

    arb_clear(beta);
}

/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_dirichlet.h"
#include "arb_hypgeom.h"

static void
_arb_pow_si(arb_t res, const arb_t x, slong y, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_set_si(a, y);
    arb_pow(res, x, a, prec);
    arb_clear(a);
}

static void
_arb_add_d(arb_t res, const arb_t x, double y, slong prec)
{
    arb_t a;
    arb_init(a);
    arb_set_d(a, y);
    arb_add(res, x, a, prec);
    arb_clear(a);
}

static void
_arb_sub_d(arb_t res, const arb_t x, double y, slong prec)
{
    _arb_add_d(res, x, -y, prec);
}

/* res = 2^((6*k+5-sigma)/4) * pi^k * (sigma + 1/2)^k * h */
static void
_pre_c_Xa(arb_t res, slong sigma, const arb_t h, ulong k, slong prec)
{
    arb_t pi, two, x1, x2;

    arb_init(pi);
    arb_init(two);
    arb_init(x1);
    arb_init(x2);

    arb_const_pi(pi, prec);
    arb_set_si(two, 2);

    arb_set_si(x1, (6*k + 5 - sigma));
    arb_mul_2exp_si(x1, x1, -2);
    arb_pow(x1, two, x1, prec);

    arb_set_si(x2, sigma);
    _arb_add_d(x2, x2, 0.5, prec);
    arb_mul(x2, x2, pi, prec);
    arb_pow_ui(x2, x2, k, prec);

    arb_mul(res, x1, x2, prec);
    arb_mul(res, res, h, prec);

    arb_clear(pi);
    arb_clear(two);
    arb_clear(x1);
    arb_clear(x2);
}

/* res = 2^((6*k+7-sigma)/4) * pi^(k - 1/2) */
static void
_pre_c_Xb(arb_t res, slong sigma, ulong k, slong prec)
{
    arb_t pi, two, x1, x2;

    arb_init(pi);
    arb_init(two);
    arb_init(x1);
    arb_init(x2);

    arb_const_pi(pi, prec);
    arb_set_si(two, 2);

    arb_set_si(x1, (6*k + 7 - sigma));
    arb_mul_2exp_si(x1, x1, -2);
    arb_pow(x1, two, x1, prec);

    arb_set_ui(x2, k);
    _arb_sub_d(x2, x2, 0.5, prec);
    arb_pow(x2, pi, x2, prec);

    arb_mul(res, x1, x2, prec);

    arb_clear(pi);
    arb_clear(two);
    arb_clear(x1);
    arb_clear(x2);
}

/*
 * res[(sigma-1)/2 - l] = binom((sigma-1)/2, l) * 2^((k+l-1)/2)
 *                        * h^(k+l+1) * gamma((k+l+1)/2, (sigma+1/2)^2 / 2*h^2)
 */
static void
_pre_c_p(arb_ptr res, slong sigma, const arb_t h, ulong k, slong prec)
{
    slong l;
    slong len = (sigma - 1)/2 + 1;
    arb_t two, f, f1, f2, u, base, x;

    arb_init(two);
    arb_init(f);
    arb_init(f1);
    arb_init(f2);
    arb_init(u);
    arb_init(base);
    arb_init(x);

    /* precompute 2^((k-1)/2) * h^(k+1) */
    arb_set_ui(two, 2);
    arb_set_si(f1, k-1);
    arb_mul_2exp_si(f1, f1, -1);
    arb_pow(f1, two, f1, prec);
    _arb_pow_si(f2, h, k+1, prec);
    arb_mul(f, f1, f2, prec);

    /* precompute (sigma + 1/2)^2 / 2*h^2 */
    arb_set_si(u, sigma);
    _arb_add_d(u, u, 0.5, prec);
    arb_div(u, u, h, prec);
    arb_sqr(u, u, prec);
    arb_mul_2exp_si(u, u, -1);

    /* precompute powers of sqrt(2)*h */
    arb_sqrt_ui(base, 2, prec);
    arb_mul(base, base, h, prec);
    _arb_vec_set_powers(res, base, len, prec);

    for (l = 0; l < len; l++)
    {
        arb_mul(res + l, res + l, f, prec);

        arb_bin_uiui(x, (sigma - 1)/2, l, prec);
        arb_mul(res + l, res + l, x, prec);

        arb_set_si(x, k + l + 1);
        arb_mul_2exp_si(x, x, -1);
        arb_hypgeom_gamma_upper(x, x, u, 0, prec);
        arb_mul(res + l, res + l, x, prec);
    }

    _arb_poly_reverse(res, res, len, len);

    arb_clear(two);
    arb_clear(f);
    arb_clear(f1);
    arb_clear(f2);
    arb_clear(u);
    arb_clear(base);
    arb_clear(x);
}

void
acb_dirichlet_platt_c_precomp_init(acb_dirichlet_platt_c_precomp_t pre,
        slong sigma, const arb_t h, ulong k, slong prec)
{
    if (!arb_is_positive(h))
    {
        flint_throw(FLINT_ERROR, "requires positive h\n");
    }
    else if (sigma % 2 == 0 || sigma < 3)
    {
        flint_throw(FLINT_ERROR, "requires odd integer sigma >= 3 (sigma=%wd)\n", sigma);
    }
    else
    {
        pre->len = (sigma - 1) / 2 + 1;
        arb_init(&pre->Xa);
        arb_init(&pre->Xb);
        pre->p = _arb_vec_init(pre->len);
        _pre_c_Xa(&pre->Xa, sigma, h, k, prec);
        _pre_c_Xb(&pre->Xb, sigma, k, prec);
        _pre_c_p(pre->p, sigma, h, k, prec);
    }
}

void
acb_dirichlet_platt_c_precomp_clear(acb_dirichlet_platt_c_precomp_t pre)
{
    arb_clear(&pre->Xa);
    arb_clear(&pre->Xb);
    _arb_vec_clear(pre->p, pre->len);
}

void
acb_dirichlet_platt_c_bound_precomp(arb_t res,
        const acb_dirichlet_platt_c_precomp_t pre, slong sigma, const arb_t t0,
        const arb_t h, slong k, slong prec)
{
    /* requires sigma + 1/2 <= t0 */
    arb_t lhs;
    arb_init(lhs);
    arb_set_si(lhs, sigma);
    _arb_add_d(lhs, lhs, 0.5, prec);
    if (!arb_le(lhs, t0))
    {
        arb_zero_pm_inf(res);
    }
    else
    {
        arb_t u, v;

        arb_init(u);
        arb_init(v);

        /* u = exp((1 + sqrt(8))/(6*t0)) */
        arb_sqrt_ui(u, 8, prec);
        arb_add_ui(u, u, 1, prec);
        arb_div_ui(u, u, 6, prec);
        arb_div(u, u, t0, prec);
        arb_exp(u, u, prec);

        /* v = (sigma + 1/2 + t0)^((sigma - 1)/2) */
        arb_add_si(v, t0, sigma, prec);
        _arb_add_d(v, v, 0.5, prec);
        _arb_pow_si(v, v, (sigma - 1)/2, prec);

        /* res = u*(v*Xa + Xb*X(t0)) */
        _arb_poly_evaluate(res, pre->p, pre->len, t0, prec);
        arb_mul(res, res, &pre->Xb, prec);
        arb_addmul(res, v, &pre->Xa, prec);
        arb_mul(res, res, u, prec);

        arb_clear(u);
        arb_clear(v);
    }
    arb_clear(lhs);
}

void
acb_dirichlet_platt_c_bound(arb_t res,
        slong sigma, const arb_t t0, const arb_t h, slong k, slong prec)
{
    acb_dirichlet_platt_c_precomp_t pre;
    acb_dirichlet_platt_c_precomp_init(pre, sigma, h, k, prec);
    acb_dirichlet_platt_c_bound_precomp(res, pre, sigma, t0, h, k, prec);
    acb_dirichlet_platt_c_precomp_clear(pre);
}

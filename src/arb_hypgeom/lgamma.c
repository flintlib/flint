/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n, const arb_t x, int use_reflect, int digamma, slong prec);
int arb_hypgeom_gamma_exact(arb_t res, const arb_t x, int reciprocal, slong prec);
void arb_hypgeom_gamma_stirling_inner(arb_t s, const arb_t z, slong N, slong prec);

void
arb_hypgeom_lgamma_stirling(arb_t y, const arb_t x, slong prec)
{
    int reflect;
    slong r, n, wp;
    arb_t t, u;
    double acc;

    /* todo: for large x (if exact or accurate enough), increase precision */
    acc = arb_rel_accuracy_bits(x);
    acc = FLINT_MAX(acc, 0);
    wp = FLINT_MIN(prec, acc + 20);
    wp = FLINT_MAX(wp, 2);
    wp = wp + FLINT_BIT_COUNT(wp);

    arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, x, 0, 0, wp);

    arb_init(t);
    arb_init(u);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    arb_init(t);
    arb_init(u);

    arb_add_ui(t, x, r, wp);
    arb_hypgeom_gamma_stirling_inner(u, t, n, wp);
    arb_hypgeom_rising_ui_rec(t, x, r, wp);
    arb_log(t, t, wp);
    arb_sub(y, u, t, prec);

    arb_clear(t);
    arb_clear(u);
}

void
arb_hypgeom_lgamma(arb_t res, const arb_t x, slong prec)
{
    if (!arb_is_positive(x) || !arb_is_finite(x))
    {
        arb_indeterminate(res);
        return;
    }

    if (arb_hypgeom_gamma_exact(res, x, 0, prec))
    {
        arb_log(res, res, prec);
        return;
    }

    if (arb_hypgeom_gamma_taylor(res, x, 0, prec))
    {
        arb_log(res, res, prec);
        return;
    }

    arb_hypgeom_lgamma_stirling(res, x, prec);
}


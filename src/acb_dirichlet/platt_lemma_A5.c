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

/* Lemma A5 requires B > h*sqrt(k) */
static int
_platt_lemma_A5_constraint(slong B, const arb_t h, slong k, slong prec)
{
    int result;
    arb_t lhs, rhs;
    arb_init(rhs);
    arb_init(lhs);
    arb_set_si(lhs, B);
    arb_sqrt_ui(rhs, (ulong) k, prec);
    arb_mul(rhs, rhs, h, prec);
    result = arb_gt(lhs, rhs);
    arb_clear(rhs);
    arb_clear(lhs);
    return result;
}


void
acb_dirichlet_platt_lemma_A5(arb_t out, slong B, const arb_t h,
        slong k, slong prec)
{
    arb_t a, b;
    arb_t x1, x2, x3, x4, x5, x6;

    if (!_platt_lemma_A5_constraint(B, h, k, prec))
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(a);
    arb_init(b);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);
    arb_init(x6);

    arb_const_pi(x1, prec);
    arb_mul_si(x1, x1, B, prec);
    arb_pow_ui(x1, x1, (ulong) k, prec);
    arb_mul_2exp_si(x1, x1, 3);

    arb_set_si(a, B);
    arb_div(a, a, h, prec);
    arb_sqr(a, a, prec);
    arb_mul_2exp_si(a, a, -3);

    arb_neg(x2, a);
    arb_exp(x2, x2, prec);

    arb_set_si(x3, 3*k - 1);
    arb_mul_2exp_si(x3, x3, -1);

    arb_set_ui(x4, 2);
    arb_pow(x4, x4, x3, prec);

    arb_set_si(b, k+1);

    arb_div_si(x5, h, B, prec);
    arb_pow_ui(x5, x5, (ulong) (k+1), prec);

    arb_mul_2exp_si(x6, b, -1);

    arb_hypgeom_gamma_upper(x6, x6, a, 0, prec);

    arb_mul(out, x4, x5, prec);
    arb_mul(out, out, x6, prec);
    arb_add(out, out, x2, prec);
    arb_mul(out, out, x1, prec);

    arb_clear(a);
    arb_clear(b);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
    arb_clear(x6);
}

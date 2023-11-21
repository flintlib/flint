/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_spherical_y(acb_t res, slong n, slong m,
    const acb_t theta, const acb_t phi, slong prec)
{
    acb_t t, u;

    if (n < 0)
    {
        if (m <= n)
        {
            acb_zero(res);
            return;
        }

        n = -1-n;
    }

    if (m > n || m < -n)
    {
        acb_zero(res);
        return;
    }

    if (n > WORD_MAX / 4)
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(t);
    acb_init(u);

    acb_sin_cos(t, u, theta, prec);

    /* P_n^m(cos(theta)) */
    acb_hypgeom_legendre_p_uiui_rec(u, n, FLINT_ABS(m), u, prec);
    acb_pow_ui(t, t, FLINT_ABS(m), prec);
    acb_mul(t, t, u, prec);

    /* exp(i m phi) */
    acb_mul_onei(u, phi);
    acb_mul_si(u, u, m, prec);
    acb_exp(u, u, prec);
    if (m < 0 && m % 2)
        acb_neg(u, u);
    acb_mul(t, t, u, prec);

    /* sqrt((2n+1)/(4pi) (n-m)!/(n+m)!) */
    arb_fac_ui(acb_realref(u), n - FLINT_ABS(m), prec);
    arb_fac_ui(acb_imagref(u), n + FLINT_ABS(m), prec);
    arb_mul_ui(acb_realref(u), acb_realref(u), 2 * n + 1, prec);
    arb_div(acb_realref(u), acb_realref(u), acb_imagref(u), prec);

    arb_const_pi(acb_imagref(u), prec);
    arb_div(acb_realref(u), acb_realref(u), acb_imagref(u), prec);
    arb_mul_2exp_si(acb_realref(u), acb_realref(u), -2);
    arb_sqrt(acb_realref(u), acb_realref(u), prec);

    acb_mul_arb(t, t, acb_realref(u), prec);

    acb_set(res, t);

    acb_clear(t);
    acb_clear(u);
}


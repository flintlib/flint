/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_calc.h"

void
acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param,
    const acb_t x, const arb_t radius, slong maxdepth, slong prec)
{
    slong i, n, depth, wp;

    arb_t pi, theta, v, s1, c1, s2, c2, st, ct;
    acb_t t, u;
    arb_t b;

    arb_init(pi);
    arb_init(theta);
    arb_init(v);

    arb_init(s1);
    arb_init(c1);
    arb_init(s2);
    arb_init(c2);
    arb_init(st);
    arb_init(ct);

    acb_init(t);
    acb_init(u);
    arb_init(b);

    wp = prec + 20;

    arb_const_pi(pi, wp);
    arb_zero_pm_inf(b);

    for (depth = 0, n = 16; depth < maxdepth; n *= 2, depth++)
    {
        arb_zero(b);

        /* theta = 2 pi / n */
        arb_div_ui(theta, pi, n, wp);
        arb_mul_2exp_si(theta, theta, 1);

        /* sine and cosine of i*theta and (i+1)*theta */
        arb_zero(s1);
        arb_one(c1);
        arb_sin_cos(st, ct, theta, wp);
        arb_set(s2, st);
        arb_set(c2, ct);

        for (i = 0; i < n; i++)
        {
            /* sine and cosine of 2 pi ([i,i+1]/n) */

            /* since we use power of two subdivision points, the
               sine and cosine are monotone on each subinterval */
            arb_union(acb_realref(t), c1, c2, wp);
            arb_union(acb_imagref(t), s1, s2, wp);
            acb_mul_arb(t, t, radius, wp);
            acb_add(t, t, x, prec);

            /* next angle */
            arb_mul(v, c2, ct, wp);
            arb_mul(c1, s2, st, wp);
            arb_sub(c1, v, c1, wp);
            arb_mul(v, c2, st, wp);
            arb_mul(s1, s2, ct, wp);
            arb_add(s1, v, s1, wp);
            arb_swap(c1, c2);
            arb_swap(s1, s2);

            func(u, t, param, 1, prec);
            acb_abs(v, u, prec);
            arb_add(b, b, v, prec);
        }

        arb_div_ui(b, b, n, prec);

        if (arb_is_positive(b))
            break;
    }

    arb_set(bound, b);

    arb_clear(pi);
    arb_clear(theta);
    arb_clear(v);

    acb_clear(t);
    acb_clear(u);
    arb_clear(b);

    arb_clear(s1);
    arb_clear(c1);
    arb_clear(s2);
    arb_clear(c2);
    arb_clear(st);
    arb_clear(ct);
}


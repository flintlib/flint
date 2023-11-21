/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sqrt(acb_t y, const acb_t x, slong prec)
{
    arb_t r, t, u;
    slong wp;
    int done;

#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    if (arb_is_zero(b))
    {
        if (arb_is_nonnegative(a))
        {
            arb_sqrt(c, a, prec);
            arb_zero(d);
            return;
        }
        else if (arb_is_nonpositive(a))
        {
            arb_neg(d, a);
            arb_sqrt(d, d, prec);
            arb_zero(c);
            return;
        }
    }

    if (arb_is_zero(a))
    {
        if (arb_is_nonnegative(b))
        {
            arb_mul_2exp_si(c, b, -1);
            arb_sqrt(c, c, prec);
            arb_set(d, c);
            return;
        }
        else if (arb_is_nonpositive(b))
        {
            arb_mul_2exp_si(c, b, -1);
            arb_neg(c, c);
            arb_sqrt(c, c, prec);
            arb_neg(d, c);
            return;
        }
    }

    wp = prec + 4;

    arb_init(r);
    arb_init(t);
    arb_init(u);

    /* r = |a+bi| */
    acb_abs(r, x, wp);

    done = 0;

    if (arf_sgn(arb_midref(a)) >= 0)
    {
        arb_add(t, r, a, wp);

        if (arb_rel_accuracy_bits(t) > 8)
        {
            /* sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i */
            arb_mul_2exp_si(u, t, 1);
            arb_sqrt(u, u, wp);
            arb_div(d, b, u, prec);

            arb_set_round(c, u, prec);
            arb_mul_2exp_si(c, c, -1);
            done = 1;
        }
        else
        {
            arb_sub(u, r, a, wp);
        }
    }
    else if (!arb_contains_zero(b))
    {
        arb_sub(u, r, a, wp);

        if (arb_rel_accuracy_bits(u) > 8)
        {
            /* sqrt(a+bi) = |b|/sqrt(2*(r-a)) + sgn(b)*sqrt((r-a)/2)*i */
            int sgn = arf_sgn(arb_midref(b));

            arb_mul_2exp_si(t, u, 1);
            arb_sqrt(t, t, wp);
            arb_div(c, b, t, prec);
            arb_abs(c, c);

            arb_set_round(d, t, prec);
            arb_mul_2exp_si(d, d, -1);
            if (sgn < 0)
                arb_neg(d, d);
            done = 1;
        }
        else
        {
            arb_add(t, r, a, wp);
        }
    }
    else
    {
        arb_add(t, r, a, wp);
        arb_sub(u, r, a, wp);
    }

    /* t = r+a, u = r-a */

    if (!done)
    {
        /* sqrt(a+bi) = sqrt((r+a)/2) + (b/|b|)*sqrt((r-a)/2)*i
                                        (sign)                      */
        arb_mul_2exp_si(t, t, -1);
        arb_mul_2exp_si(u, u, -1);
        arb_sqrtpos(c, t, prec);

        if (arb_is_nonnegative(b))
        {
            arb_sqrtpos(d, u, prec);
        }
        else if (arb_is_nonpositive(b))
        {
            arb_sqrtpos(d, u, prec);
            arb_neg(d, d);
        }
        else
        {
            arb_sqrtpos(t, u, wp);
            arb_neg(u, t);
            arb_union(d, t, u, prec);
        }
    }

    arb_clear(r);
    arb_clear(t);
    arb_clear(u);

#undef a
#undef b
#undef c
#undef d
}

void
acb_sqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (analytic && arb_contains_zero(acb_imagref(z)) &&
        !arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_sqrt(res, z, prec);
    }
}


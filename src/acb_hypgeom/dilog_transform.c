/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_dilog_transform(acb_t res, const acb_t z, int algorithm, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    if (algorithm == 1)
    {
        /* Li_2(z) = -Li_2(1/z) - log(-z)^2/2 - pi^2/6,  z not in (0,1) */
        arf_set_ui_2exp_si(arb_midref(acb_realref(t)), 1, -1);
        mag_set_ui_2exp_si(arb_radref(acb_realref(t)), 1, -1);

        if (acb_overlaps(z, t))
        {
            acb_indeterminate(res);
        }
        else
        {
            acb_inv(t, z, prec);
            acb_hypgeom_dilog_zero(t, t, prec);
            acb_neg(u, z);
            acb_log(u, u, prec);
            acb_mul(u, u, u, prec);
            acb_mul_2exp_si(u, u, -1);
            acb_add(t, t, u, prec);
            acb_const_pi(u, prec);
            acb_mul(u, u, u, prec);
            acb_div_ui(u, u, 6, prec);
            acb_add(t, t, u, prec);
            acb_neg(res, t);
        }
    }
    else if (algorithm == 2)
    {
        /* Li_2(z) = -Li_2(1-z) - log(1-z) log(z) + pi^2/6 */
        if (acb_is_one(z))
        {
            acb_zero(res);
        }
        else
        {
            acb_sub_ui(t, z, 1, prec);
            acb_neg(t, t);
            acb_hypgeom_dilog_zero(u, t, prec);
            acb_log(t, t, prec);
            acb_log(res, z, prec);
            acb_mul(res, res, t, prec);
            acb_add(res, res, u, prec);
        }

        acb_const_pi(t, prec);
        acb_mul(t, t, t, prec);
        acb_div_ui(t, t, 6, prec);
        acb_sub(res, t, res, prec);
    }
    else if (algorithm == 3)
    {
        /* Li_2(z) = -Li_2(z/(z-1)) - log(1-z)^2/2,  z not in (1,inf) */
        acb_sub_ui(t, z, 1, prec);

        if (!arb_is_negative(acb_realref(t)))
        {
            acb_indeterminate(res);
        }
        else
        {
            acb_div(u, z, t, prec);
            acb_hypgeom_dilog_zero(u, u, prec);
            acb_neg(t, t);
            acb_log(t, t, prec);
            acb_mul(t, t, t, prec);
            acb_mul_2exp_si(t, t, -1);
            acb_add(t, t, u, prec);
            acb_neg(res, t);
        }
    }
    else if (algorithm == 4)
    {
        /* Li_2(z) = Li_2(1/(1-z)) + log(1-z) [log(1-z)/2 - log(-z)] - pi^2/6 */
        acb_sub_ui(t, z, 1, prec);
        acb_neg(t, t);

        acb_inv(u, t, prec);
        acb_hypgeom_dilog_zero(u, u, prec);

        acb_log(t, t, prec);
        acb_neg(res, z);
        acb_log(res, res, prec);
        acb_mul_2exp_si(res, res, 1);
        acb_sub(res, t, res, prec);
        acb_mul_2exp_si(res, res, -1);
        acb_addmul(u, res, t, prec);

        acb_const_pi(t, prec);
        acb_mul(t, t, t, prec);
        acb_div_ui(t, t, 6, prec);
        acb_sub(res, u, t, prec);
    }
    else if (algorithm >= 5 && algorithm <= 7)
    {
        if (arb_contains_zero(acb_imagref(z)))
        {
            acb_indeterminate(res);
        }
        else
        {
            acb_t a;
            acb_init(a);

            if (algorithm == 5)
            {
                acb_onei(a);
                /* Li_2(i) = -pi^2/48 + C i */
                arb_const_catalan(acb_imagref(u), prec);
                arb_const_pi(acb_realref(u), prec);
                arb_mul(acb_realref(u), acb_realref(u), acb_realref(u), prec);
                arb_div_si(acb_realref(u), acb_realref(u), -48, prec);
            }
            else if (algorithm == 6)
            {
                /* Li_2((1+i)/2) = (5 pi^2 / 96 - log(2)^2/8) + (C - pi log(2) / 8) i */
                arb_t t;
                arb_init(t);
                acb_set_d_d(a, 0.5, 0.5);
                arb_const_pi(t, prec);
                arb_const_log2(acb_imagref(u), prec);
                arb_mul(acb_realref(u), acb_imagref(u), acb_imagref(u), prec);
                arb_mul(acb_imagref(u), acb_imagref(u), t, prec);
                acb_mul_2exp_si(u, u, -3);
                arb_mul(t, t, t, prec);
                arb_mul_ui(t, t, 5, prec);
                arb_div_ui(t, t, 96, prec);
                arb_sub(acb_realref(u), t, acb_realref(u), prec);
                arb_const_catalan(t, prec);
                arb_sub(acb_imagref(u), t, acb_imagref(u), prec);
                arb_clear(t);
            }
            else
            {
                /* Li_2(1+i) = pi^2/16 + (C + pi log(2)/4) i */
                arb_t t;
                arb_init(t);
                acb_set_d_d(a, 1.0, 1.0);
                arb_const_pi(acb_realref(u), prec);
                arb_mul_2exp_si(acb_realref(u), acb_realref(u), -2);
                arb_const_log2(t, prec);
                arb_mul(acb_imagref(u), acb_realref(u), t, prec);
                arb_const_catalan(t, prec);
                arb_add(acb_imagref(u), acb_imagref(u), t, prec);
                arb_mul(acb_realref(u), acb_realref(u), acb_realref(u), prec);
                arb_clear(t);
            }

            if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
            {
                acb_conj(a, a);
                acb_conj(u, u);
            }

            acb_hypgeom_dilog_bitburst(res, t, z, prec);
            acb_add(res, res, u, prec);
            acb_hypgeom_dilog_continuation(t, a, t, prec);
            acb_add(res, res, t, prec);

            acb_clear(a);
        }
    }
    else
    {
        flint_throw(FLINT_ERROR, "unknown algorithm\n");
    }

    acb_clear(t);
    acb_clear(u);
}


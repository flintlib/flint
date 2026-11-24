/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "arb_hypgeom/impl.h"
#include "acb_calc.h"

/* a random dyadic number in [-1, 1] times scale */
static void
randd(arb_t x, flint_rand_t state, double scale)
{
    arb_set_si(x, (slong) n_randint(state, 2001) - 1000);
    arb_div_ui(x, x, 1000, 30);
    arb_get_mid_arb(x, x);
    arb_mul_2exp_si(x, x, (slong) floor(log2(scale) + 0.5));
}

typedef struct
{
    arb_srcptr lam, mu, al, be;
    int m;
}
kernel_t;

/* exp(lambda t + sum_j mu_j log(al_j + be_j t)) */
static int
integrand(acb_ptr out, const acb_t t, void * param, slong order, slong prec)
{
    const kernel_t * K = param;
    acb_t s, u;
    int j;

    if (order > 1)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    acb_init(s);
    acb_init(u);

    acb_mul_arb(s, t, K->lam, prec);
    for (j = 0; j < K->m; j++)
    {
        acb_mul_arb(u, t, K->be + j, prec);
        acb_add_arb(u, u, K->al + j, prec);
        /* (analytic where the real part is positive) */
        if (order == 1 && !arb_is_positive(acb_realref(u)))
            acb_indeterminate(u);
        else
            acb_log(u, u, prec);
        acb_addmul_arb(s, u, K->mu + j, prec);
    }
    acb_exp(out, s, prec);

    acb_clear(s);
    acb_clear(u);
    return 0;
}

TEST_FUNCTION_START(arb_hypgeom_integrate_exp_log, state)
{
    slong iter, computed = 0;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t lam, r1, r2, r3, t, u;
        arb_ptr mu, al, be, mu2;
        slong prec1, prec2;
        int m, j, type, ok1, ok2, have3 = 0, sing0;
        double B, scale;

        arb_init(lam); arb_init(r1); arb_init(r2); arb_init(r3); arb_init(t); arb_init(u);
        mu = _arb_vec_init(3);
        mu2 = _arb_vec_init(3);
        al = _arb_vec_init(3);
        be = _arb_vec_init(3);

        prec1 = 20 + n_randint(state, (n_randint(state, 4) == 0) ? 400 : 100);
        prec2 = prec1 + 20 + n_randint(state, 100);
        m = 1 + n_randint(state, 3);
        type = n_randint(state, 4);
        scale = exp(log(1000.0) * n_randint(state, 1001) / 1000.0);
        if (type == 0)
            scale = FLINT_MIN(scale, 30.0);

        if (type <= 1)
            B = ldexp(1.0 + n_randint(state, 50), -(int) n_randint(state, 6));
        else
            B = INFINITY;

        randd(lam, state, scale);
        if (type >= 2)
        {
            /* lambda < 0 */
            arb_abs(lam, lam);
            arb_add_ui(lam, lam, 1, 60);
            arb_mul_2exp_si(lam, lam, -(slong) n_randint(state, 4));
            arb_neg(lam, lam);
        }

        for (j = 0; j < m; j++)
        {
            randd(mu + j, state, scale);
            if (n_randint(state, 8) == 0)
                arb_zero(mu + j);
            if (n_randint(state, 4) == 0)
                arb_floor(mu + j, mu + j, 60);

            /* al_j > 0, and be_j such that al_j + be_j t > 0 */
            randd(al + j, state, 4.0);
            arb_abs(al + j, al + j);
            arb_add_ui(al + j, al + j, 1, 60);
            arb_mul_2exp_si(al + j, al + j, -1);
            randd(be + j, state, 4.0);
            if (type >= 2)
                arb_abs(be + j, be + j);
            if (type >= 2 && arb_is_zero(be + j))
                arb_one(be + j);
            if (type <= 1 && arb_is_negative(be + j))
            {
                /* al_j / B / 2 at most */
                arb_set_d(t, 2.0 * B);
                arb_div(t, al + j, t, 60);
                arb_get_mid_arb(t, t);
                arb_neg(t, t);
                arb_max(be + j, be + j, t, 60);
            }
            if (n_randint(state, 3) == 0)
            {
                arb_set_si(be + j, arb_is_negative(be + j) ? -1 : 1);
                if (type <= 1 && arb_is_negative(be + j))
                {
                    arb_set_d(t, 2.0 * B);
                    arb_add(al + j, al + j, t, 60);
                }
            }
        }

        /* a term which vanishes at 0 */
        sing0 = -1;
        if (type != 0)
        {
            sing0 = 0;
            arb_zero(al);
            arb_abs(be, be);
            if (arb_is_zero(be))
                arb_one(be);
            arb_abs(mu, mu);
        }

        /* inexact parameters */
        _arb_vec_set(mu2, mu, m);
        if (n_randint(state, 3) == 0)
        {
            for (j = 0; j < m; j++)
                if (j != sing0 || n_randint(state, 2))
                    if (!arb_is_zero(mu2 + j))
                        mag_set_ui_2exp_si(arb_radref(mu2 + j), 1, -20 - (slong) n_randint(state, 60));
        }

        ok1 = _arb_hypgeom_integrate_exp_log(r1, lam, mu2, al, be, m, sing0, B, -INFINITY, prec1);
        ok2 = _arb_hypgeom_integrate_exp_log(r2, lam, mu, al, be, m, sing0, B, -INFINITY, prec2);
        computed += ok1;

        if ((ok1 && !arb_is_finite(r1)) || (ok2 && !arb_is_finite(r2)) || (ok1 && ok2 && !arb_overlaps(r1, r2)))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("type %d, m = %d, B = %g, prec1 = %wd, prec2 = %wd\n", type, m, B, prec1, prec2);
            flint_printf("lam = "); arb_printd(lam, 20); flint_printf("\n");
            for (j = 0; j < m; j++)
            {
                flint_printf("mu = "); arb_printd(mu2 + j, 20);
                flint_printf(", al = "); arb_printd(al + j, 20);
                flint_printf(", be = "); arb_printd(be + j, 20); flint_printf("\n");
            }
            flint_printf("r1 = "); arb_printd(r1, 50); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printd(r2, 50); flint_printf("\n\n");
            flint_abort();
        }

        /* reference values */
        if (type == 0)
        {
            /* by acb_calc_integrate */
            kernel_t K;
            acb_t a, b, I;
            mag_t tol;

            acb_init(a); acb_init(b); acb_init(I);
            mag_init(tol);

            K.lam = lam; K.mu = mu; K.al = al; K.be = be; K.m = m;
            acb_set_d(b, B);
            mag_set_ui_2exp_si(tol, 1, -1000);
            acb_calc_integrate(I, integrand, &K, a, b, 60, tol, NULL, 80);
            arb_set(r3, acb_realref(I));
            have3 = arb_is_finite(r3) && arb_rel_accuracy_bits(r3) > 20;

            acb_clear(a); acb_clear(b); acb_clear(I);
            mag_clear(tol);
        }
        else if (type == 1 && m == 1)
        {
            /* int_0^B e^(lam t) (be t)^mu dt = be^mu B^(mu + 1) M(mu + 1, mu + 2, lam B) / (mu + 1) */
            arb_set_d(t, B);
            arb_mul(u, lam, t, prec2 + 30);
            arb_add_ui(r3, mu, 1, prec2 + 30);
            arb_pow(t, t, r3, prec2 + 30);
            arb_div(t, t, r3, prec2 + 30);
            arb_add_ui(r2, r3, 1, prec2 + 30);
            arb_hypgeom_1f1(r3, r3, r2, u, 0, prec2 + 30);
            arb_mul(r3, r3, t, prec2 + 30);
            arb_pow(t, be, mu, prec2 + 30);
            arb_mul(r3, r3, t, prec2 + 30);
            have3 = arb_is_finite(r3) && arb_rel_accuracy_bits(r3) > 20;
        }
        else if (type >= 2 && m == 1)
        {
            /* int_0^inf e^(lam t) (be t)^mu dt = be^mu Gamma(mu + 1) / (-lam)^(mu + 1) */
            arb_add_ui(t, mu, 1, prec2 + 30);
            arb_gamma(r3, t, prec2 + 30);
            arb_neg(u, lam);
            arb_pow(u, u, t, prec2 + 30);
            arb_div(r3, r3, u, prec2 + 30);
            arb_pow(t, be, mu, prec2 + 30);
            arb_mul(r3, r3, t, prec2 + 30);
            have3 = arb_is_finite(r3);
        }

        if (have3 && ok1 && !arb_overlaps(r1, r3))
        {
            flint_printf("FAIL: reference value\n\n");
            flint_printf("type %d, m = %d, B = %g, prec1 = %wd, prec2 = %wd\n", type, m, B, prec1, prec2);
            flint_printf("lam = "); arb_printd(lam, 20); flint_printf("\n");
            for (j = 0; j < m; j++)
            {
                flint_printf("mu = "); arb_printd(mu2 + j, 20);
                flint_printf(", al = "); arb_printd(al + j, 20);
                flint_printf(", be = "); arb_printd(be + j, 20); flint_printf("\n");
            }
            flint_printf("r1 = "); arb_printd(r1, 50); flint_printf("\n\n");
            flint_printf("r3 = "); arb_printd(r3, 50); flint_printf("\n\n");
            flint_abort();
        }

        /* the closed forms are computed */
        if (m == 1 && type >= 1 && !ok1)
        {
            flint_printf("FAIL: not computed\n\n");
            flint_printf("type %d, B = %g, prec1 = %wd\n", type, B, prec1);
            flint_printf("lam = "); arb_printd(lam, 20); flint_printf("\n");
            flint_printf("mu = "); arb_printd(mu2, 20);
            flint_printf(", be = "); arb_printd(be, 20); flint_printf("\n");
            flint_abort();
        }

        arb_clear(lam); arb_clear(r1); arb_clear(r2); arb_clear(r3); arb_clear(t); arb_clear(u);
        _arb_vec_clear(mu, 3);
        _arb_vec_clear(mu2, 3);
        _arb_vec_clear(al, 3);
        _arb_vec_clear(be, 3);
    }

    if (iter >= 20 && computed < 0.9 * iter)
    {
        flint_printf("FAIL: %wd of %wd integrals computed\n\n", computed, iter);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}

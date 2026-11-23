/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_u_asymp_proper(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_pow(t, z, a, prec);
    acb_hypgeom_u_asymp(res, a, b, z, -1, prec);
    acb_div(res, res, t, prec);
    acb_clear(t);
}

static void
_acb_add_error_2exp_si(acb_t x, slong e)
{
    arb_add_error_2exp_si(acb_realref(x), e);
    arb_add_error_2exp_si(acb_imagref(x), e);
}

TEST_FUNCTION_START(acb_hypgeom_u, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a0, a1, a2, b, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;

        acb_init(a0);
        acb_init(a1);
        acb_init(a2);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 700);
        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_param(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        if (prec0 <= 300 && prec1 <= 300 && prec2 <= 300 && n_randint(state, 2))
        {
            arb_zero(acb_imagref(a0));
            arb_zero(acb_imagref(b));
            arb_zero(acb_imagref(z));
        }

        acb_add_ui(a1, a0, 1, prec0);
        acb_add_ui(a2, a0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w0, a0, b, z, prec0);
                break;
            case 1:
                acb_hypgeom_u_1f1(w0, a0, b, z, prec0);
                break;
            default:
                acb_hypgeom_u(w0, a0, b, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w1, a0, b, z, prec1);
                break;
            case 1:
                acb_hypgeom_u_1f1(w1, a0, b, z, prec1);
                break;
            default:
                acb_hypgeom_u(w1, a0, b, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w1, a1, b, z, prec1);
                break;
            case 1:
                acb_hypgeom_u_1f1(w1, a1, b, z, prec1);
                break;
            default:
                acb_hypgeom_u(w1, a1, b, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w2, a2, b, z, prec2);
                break;
            case 1:
                acb_hypgeom_u_1f1(w2, a2, b, z, prec2);
                break;
            default:
                acb_hypgeom_u(w2, a2, b, z, prec2);
        }

        acb_set(t, w0);

        acb_mul_2exp_si(u, a0, 1);
        acb_sub(u, u, b, prec0);
        acb_add(u, u, z, prec0);
        acb_add_ui(u, u, 2, prec0);

        acb_submul(t, w1, u, prec0);

        acb_sub(u, a2, b, prec0);
        acb_mul(u, u, a1, prec0);
        acb_addmul(t, w2, u, prec0);

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(a2);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    /* Polynomials in 1/z (a or a - b + 1 a nonpositive integer) of
       large degree whose terms cancel: accurate at the precision of a
       double (the integrals, or more precision for the polynomial). */
    {
        static const double par[5][3] = { { -1000.0, 1.0, 10.0 }, { 1000.0, 2000.0, 1000.0 },
            { 10000.0, 20000.0, 10000.0 }, { 10.0, 500.0, 100.0 }, { -200.5, -100.5, 30.0 } };
        acb_t a, b, z, w1, w2;
        slong i;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        /* (from 5: inexact input which is much more accurate than the
           precision; 10 and 11 are not polynomials) */
        for (i = 0; i < 12; i++)
        {
            if (i < 10)
            {
                acb_set_d(a, par[i % 5][0]);
                acb_set_d(b, par[i % 5][1]);
                acb_set_d(z, par[i % 5][2]);
            }
            else
            {
                acb_set_d(a, (i == 10) ? -300.5 : 30.25);
                acb_set_d(b, 0.375);
                acb_set_d(z, (i == 10) ? 50.0 : 5.0);
            }

            acb_hypgeom_u(w2, a, b, z, 4000);

            if (i >= 5)
            {
                mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -3000);

                if (i % 2)
                {
                    mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -3000);
                    mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -3000);
                }
            }

            acb_hypgeom_u(w1, a, b, z, 64);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 3000
                || acb_rel_accuracy_bits(w1) < 50)
            {
                flint_printf("FAIL: polynomial, case %wd\n\n", i);
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    /* b close to an integer n (the expansion in b at n): points, balls
       which contain n and balls which do not, against the connection
       formula at the midpoints with the precision that its cancellation
       needs, and against the limit at n. */
    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, z, am, bm, zm, w, r, s, t, u;
        slong prec, wp, n, e;
        int contains, mode;

        acb_init(a); acb_init(b); acb_init(z);
        acb_init(am); acb_init(bm); acb_init(zm);
        acb_init(w); acb_init(r); acb_init(s); acb_init(t); acb_init(u);

        prec = 30 + n_randint(state, 300);
        n = (slong) n_randint(state, 15) - 7;

        /* 0: a point (the expansion for about half of the distances
           2^-e); 1: a ball without n whose radius is much smaller than
           the error bound of the expansion; 2: a ball with n; 3: any
           of these by the main function */
        mode = n_randint(state, 4);
        if (mode == 0)
            e = 13 + n_randint(state, 2 * prec + 30);
        else if (mode == 1)
            e = 13 + n_randint(state, prec + 30);
        else
            e = 13 + n_randint(state, 150);

        /* a: also nonpositive integers and integers a - n + 1 <= 0 */
        switch (n_randint(state, 4))
        {
            case 0:
                acb_set_si(a, (slong) n_randint(state, 12) - 8);
                break;
            case 1:
                arb_set_si(acb_realref(a), (slong) n_randint(state, 400) - 200);
                arb_mul_2exp_si(acb_realref(a), acb_realref(a), -4);
                break;
            default:
                arb_set_si(acb_realref(a), (slong) n_randint(state, 400) - 200);
                arb_set_si(acb_imagref(a), (slong) n_randint(state, 400) - 200);
                acb_mul_2exp_si(a, a, -4);
        }

        arb_set_ui(acb_realref(z), 1 + n_randint(state, 300));
        if (n_randint(state, 2))
            arb_set_si(acb_imagref(z), (slong) n_randint(state, 400) - 200);
        if (n_randint(state, 4) == 0)
            arb_neg(acb_realref(z), acb_realref(z));
        acb_mul_2exp_si(z, z, -4);

        acb_set_si(b, n);
        switch (n_randint(state, 4))
        {
            case 0:
                arb_one(acb_realref(t));
                arb_one(acb_imagref(t));
                break;
            case 1:
                acb_set_si(t, -1);
                break;
            default:
                acb_one(t);
        }
        acb_mul_2exp_si(t, t, -e);
        acb_add(b, b, t, 1000);

        acb_set(am, a);
        acb_set(bm, b);
        acb_set(zm, z);

        contains = 0;
        switch ((mode == 3) ? n_randint(state, 4) : mode)
        {
            case 1:
                _acb_add_error_2exp_si(b, -e - prec - 30);
                break;
            case 2:
                contains = 1;
                _acb_add_error_2exp_si(b, -e + 2);
                break;
            case 3:
                _acb_add_error_2exp_si(b, -e - 5 - (slong) n_randint(state, 100));
                break;
            default:
                break;
        }

        if (mode >= 2 && n_randint(state, 4) == 0)
            _acb_add_error_2exp_si(a, -20 - (slong) n_randint(state, 100));
        if (mode >= 2 && n_randint(state, 4) == 0)
            _acb_add_error_2exp_si(z, -20 - (slong) n_randint(state, 100));

        if (mode == 3)
            acb_hypgeom_u(w, a, b, z, prec);
        else
            acb_hypgeom_u_1f1(w, a, b, z, prec);

        /* pi / sin(pi b) (1F~1(a, b, z) / Gamma(a - b + 1)
            - z^(1-b) 1F~1(a - b + 1, 2 - b, z) / Gamma(a)) */
        wp = prec + 2 * e + 100;
        acb_sub(t, am, bm, wp);
        acb_add_ui(t, t, 1, wp);
        acb_rgamma(u, t, wp);
        acb_hypgeom_pfq(r, am, 1, bm, 1, zm, 1, wp);
        acb_mul(r, r, u, wp);
        acb_sub_ui(u, bm, 2, wp);
        acb_neg(u, u);
        acb_hypgeom_pfq(s, t, 1, u, 1, zm, 1, wp);
        acb_sub_ui(u, bm, 1, wp);
        acb_neg(u, u);
        acb_pow(u, zm, u, wp);
        acb_mul(s, s, u, wp);
        acb_rgamma(u, am, wp);
        acb_mul(s, s, u, wp);
        acb_sub(r, r, s, wp);
        acb_sin_pi(u, bm, wp);
        acb_div(r, r, u, wp);
        acb_const_pi(u, wp);
        acb_mul(r, r, u, wp);

        if (!acb_is_finite(r) || !acb_overlaps(w, r))
        {
            flint_printf("FAIL: b close to an integer\n\n");
            flint_printf("prec = %wd, e = %wd\n\n", prec, e);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 60); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (contains)
        {
            acb_set_si(t, n);
            acb_hypgeom_u_1f1(s, am, t, zm, prec + 30);

            if (!acb_overlaps(w, s))
            {
                flint_printf("FAIL: b close to an integer, the limit\n\n");
                flint_printf("prec = %wd, e = %wd\n\n", prec, e);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 60); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(z);
        acb_clear(am); acb_clear(bm); acb_clear(zm);
        acb_clear(w); acb_clear(r); acb_clear(s); acb_clear(t); acb_clear(u);
    }

    /* Accuracy for an integer b and for b close to an integer: the
       two series of the limit by their estimated size (one of them has
       a lower parameter which is a nonpositive integer), the
       asymptotic series with more precision at a turning point, and
       the expansion in b. */
    {
        /* a, b - n, z; b = n + 2^e (e = 0: b = n; rad: b = n +/- 2^-33) */
        static const double par[9][4] = { { 5.0, 0.0, 2.0, 3.0 }, { 1.0, -10.0, 2.0, 0.0 },
            { -1000.5, 0.0, 3.0, 4008.0 }, { 1.0 / 3, 0.0, 1.0, 1.0 }, { 1.0 / 3, 0.0, -3.0, 2.0 },
            { 0.25, 0.0, 1.0, 1.0 }, { 0.25, 0.0, 2.0, 0.5 }, { 0.25, 0.0, 1.0, 1.0 },
            { 0.25, 0.0, 1.0, 1.0 } };
        static const slong expo[9] = { 0, 0, 0, 0, 0, -30, -50, -332, 0 };
        acb_t a, b, z, w1, w2;
        slong i, prec;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 9; i++)
        {
            acb_set_d_d(a, par[i][0], par[i][1]);
            acb_set_d(b, par[i][2]);
            acb_set_d(z, par[i][3]);

            if (i == 1)
            {
                /* U(1 - 10i, 2, 40i) */
                acb_mul_onei(z, z);
                arb_set_d(acb_imagref(z), 40.0);
            }

            if (expo[i] != 0)
            {
                acb_one(w1);
                acb_mul_2exp_si(w1, w1, expo[i]);
                acb_add(b, b, w1, 1000);
            }

            acb_hypgeom_u(w2, a, b, z, 1500);

            for (prec = 64; prec <= 333; prec += 269)
            {
                slong goal = prec - 12;

                if (i == 8)
                {
                    /* b = 1 +/- 2^-33: the derivative in b is of the
                       order of 1 */
                    acb_one(b);
                    _acb_add_error_2exp_si(b, -33);
                    goal = 28;
                }

                acb_hypgeom_u(w1, a, b, z, prec);

                if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 1400
                    || acb_rel_accuracy_bits(w1) < goal)
                {
                    flint_printf("FAIL: integer b, case %wd, prec = %wd\n\n", i, prec);
                    flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                    flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    /* A parameter rounded at the precision, where the series cancel by
       more bits than the input has: the two integrals on the input, if
       they cost at most four times the series at the precision at
       which they succeed. */
    {
        static const double par[3][3] = { { -100.0, 0.5, 50.0 }, { -10000.0, 1.0, 10000.0 },
            { -1000.0, 2000.5, 1000.0 } };
        acb_t a, b, z, w1, w2;
        slong i;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 3; i++)
        {
            /* a = n - 1/3 */
            acb_set_ui(a, 1);
            acb_div_ui(a, a, 3, 64);
            acb_neg(a, a);
            acb_add_si(a, a, (slong) par[i][0], 64);
            acb_set_d(b, par[i][1]);
            acb_set_d(z, par[i][2]);

            acb_hypgeom_u(w1, a, b, z, 64);

            acb_get_mid(a, a);
            acb_hypgeom_u(w2, a, b, z, 300);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 200
                || acb_rel_accuracy_bits(w1) < 36)
            {
                flint_printf("FAIL: rounded input, case %wd\n\n", i);
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_m, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a0, a1, a2, b, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;
        int regularized, ebits;

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

        if (n_randint(state, 5) == 0)
            ebits = 100;
        else
            ebits = 10;

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest_param(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        regularized = n_randint(state, 2);

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
                acb_hypgeom_m_asymp(w0, a0, b, z, regularized, prec0);
                break;
            case 1:
                acb_hypgeom_m_1f1(w0, a0, b, z, regularized, prec0);
                break;
            default:
                acb_hypgeom_m(w0, a0, b, z, regularized, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w1, a0, b, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a0, b, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a0, b, z, regularized, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("regularized = %d\n\n", regularized);
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
                acb_hypgeom_m_asymp(w1, a1, b, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a1, b, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a1, b, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w2, a2, b, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_m_1f1(w2, a2, b, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_m(w2, a2, b, z, regularized, prec2);
        }

        /* 1F1(a,b,z) * (a-b+1) */
        acb_sub(u, a1, b, prec0);
        acb_mul(t, w0, u, prec0);

        /* 1F1(a+1,b,z) * (2a - b + z + 2) */
        acb_mul_2exp_si(u, a0, 1);
        acb_sub(u, u, b, prec0);
        acb_add(u, u, z, prec0);
        acb_add_ui(u, u, 2, prec0);
        acb_submul(t, w1, u, prec0);

        /* 1F1(a+2,b,z) * -(a+1) */
        acb_addmul(t, w2, a1, prec0);

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

        acb_add_ui(a1, b, 1, prec0);
        acb_add_ui(a2, b, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w1, a0, a1, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a0, a1, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a0, a1, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w2, a0, a2, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_m_1f1(w2, a0, a2, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_m(w2, a0, a2, z, regularized, prec2);
        }

        if (regularized)
        {
            acb_set(t, w0);

            acb_add(u, b, z, prec0);
            acb_submul(t, w1, u, prec0);

            acb_sub(u, a0, a1, prec0);
            acb_mul(u, u, z, prec0);
            acb_submul(t, w2, u, prec0);
        }
        else
        {
            acb_mul(t, w0, b, prec0);
            acb_mul(t, t, a1, prec0);

            acb_add(u, b, z, prec0);
            acb_mul(u, u, a1, prec0);
            acb_submul(t, w1, u, prec0);

            acb_sub(u, a0, a1, prec0);
            acb_mul(u, u, z, prec0);
            acb_submul(t, w2, u, prec0);
        }

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation 2\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m(a0, a0, b, z, regularized, prec0);
                acb_set(t, a0);
                break;
            case 1:
                acb_hypgeom_m(b, a0, b, z, regularized, prec0);
                acb_set(t, b);
                break;
            default:
                acb_hypgeom_m(z, a0, b, z, regularized, prec2);
                acb_set(t, z);
        }

        if (!acb_overlaps(t, w0))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
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

    /* Inexact input which is much more accurate than the precision,
       where the series cancel: as accurate as for exact input. */
    {
        static const double par[4][3] = { { 166.625, 0.375, -300.0 }, { -300.0, 0.375, 50.0 },
            { -300.5, 0.375, 50.0 }, { 40.25, 80.5, -200.0 } };
        acb_t a, b, z, w1, w2;
        slong i, p;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 8; i++)
        {
            p = (i < 4) ? 64 : 333;

            acb_set_d(a, par[i % 4][0]);
            acb_set_d(b, par[i % 4][1]);
            acb_set_d(z, par[i % 4][2]);

            acb_hypgeom_m(w2, a, b, z, 0, 3000);

            mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -2000);
            mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -2000);
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -2000);

            acb_hypgeom_m(w1, a, b, z, 0, p);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 2000
                || acb_rel_accuracy_bits(w1) < p - 14)
            {
                flint_printf("FAIL: accurate input, case %wd\n\n", i);
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

    /* A parameter rounded at the precision, where the series cancel by
       more bits than the input has: the integrals on the input, if
       they cost at most four times the series at the precision at
       which they succeed (the two integrals for U; 3 is the case in
       which the value is computed by them at this precision). */
    {
        static const double par[4][4] = { { -1000.0, 2000.5, 1000.0, 64 },
            { -1000.0, 2000.5, 535.0, 64 }, { -10000.0, 1.0, 10000.0, 64 },
            { -10000.0, 1.0, 10000.0, 333 } };
        acb_t a, b, z, w1, w2;
        slong i, p;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 4; i++)
        {
            p = (slong) par[i][3];

            /* a = n - 1/3 */
            acb_set_ui(a, 1);
            acb_div_ui(a, a, 3, p);
            acb_neg(a, a);
            acb_add_si(a, a, (slong) par[i][0], p);
            acb_set_d(b, par[i][1]);
            acb_set_d(z, par[i][2]);

            acb_hypgeom_m(w1, a, b, z, 0, p);

            acb_get_mid(a, a);
            acb_hypgeom_m(w2, a, b, z, 0, 3 * p + 100);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 2 * p
                || acb_rel_accuracy_bits(w1) < p - 28)
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

    /* Integers 0 < a < b: the finite asymptotic expansion where it is
       shorter than the series, with the bits which its two parts lose
       for small |z|; against the series. */
    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, z, w1, w2;
        slong prec, e;

        acb_init(a);
        acb_init(b);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        prec = 30 + n_randint(state, 400);

        acb_set_ui(a, 1 + n_randint(state, 40));
        acb_add_ui(b, a, 1 + n_randint(state, 40), prec);

        /* |z| from 2^-40 to 2^12 */
        e = (slong) n_randint(state, 53) - 40;
        arb_set_si(acb_realref(z), (slong) n_randint(state, 2000) - 1000);
        if (n_randint(state, 2))
            arb_set_si(acb_imagref(z), (slong) n_randint(state, 2000) - 1000);
        if (acb_is_zero(z))
            acb_one(z);
        acb_mul_2exp_si(z, z, e - 10);

        acb_hypgeom_m_1f1(w2, a, b, z, 0, 2 * prec + 6000);

        if (n_randint(state, 2) == 0)
        {
            acb_hypgeom_m(w1, a, b, z, 0, prec);
        }
        else
        {
            acb_t g;
            acb_init(g);
            acb_hypgeom_m(w1, a, b, z, 1, prec);
            acb_rgamma(g, b, 2 * prec + 6000);
            acb_mul(w2, w2, g, 2 * prec + 6000);
            acb_clear(g);
        }

        if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < prec + 20
            || acb_rel_accuracy_bits(w1) < prec - 12)
        {
            flint_printf("FAIL: integers 0 < a < b\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    TEST_FUNCTION_END(state);
}

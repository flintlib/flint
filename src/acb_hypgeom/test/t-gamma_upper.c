/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

int _mag_gt_norm_ui(const mag_t a, const mag_t b, const mag_t c, ulong n);

static void
_accuracy_regression_test(const acb_t s, const acb_t z,
        int regularized, slong prec, slong issue, slong accuracy)
{
    acb_t g;
    acb_init(g);
    acb_hypgeom_gamma_upper(g, s, z, regularized, prec);
    if (acb_rel_accuracy_bits(g) < accuracy)
    {
        flint_printf("FAIL: accuracy regression in issue #%wd\n\n", issue);
        flint_printf("prec = %wd\n\n", prec);
        flint_printf("regularized = %d\n\n", regularized);
        flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
        flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
        flint_printf("g = "); acb_printd(g, 30); flint_printf("\n\n");
        flint_abort();
    }
    acb_clear(g);
}

TEST_FUNCTION_START(acb_hypgeom_gamma_upper, state)
{
    slong iter;

    /* special accuracy test -- see nemo #38 */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, z, res;
        slong prec, goal;
        int regularized;

        acb_init(a);
        acb_init(z);
        acb_init(res);

        acb_set_si(a, n_randint(state, 100) - 50);

        do {
            acb_set_si(z, n_randint(state, 100) - 50);
        } while (acb_is_zero(z));

        regularized = n_randint(state, 3);

        goal = 2 + n_randint(state, 4000);

        for (prec = 2 + n_randint(state, 1000); ; prec *= 2)
        {
            acb_hypgeom_gamma_upper(res, a, z, regularized, prec);

            if (acb_rel_accuracy_bits(res) > goal)
                break;

            if (prec > 10000)
            {
                printf("FAIL (convergence)\n");
                flint_printf("regularized = %d\n\n", regularized);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("res = "); acb_printd(res, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(z);
        acb_clear(res);
    }

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a0, a1, b, z, w0, w1, t, u;
        slong prec0, prec1;
        int regularized;

        acb_init(a0);
        acb_init(a1);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(t);
        acb_init(u);

        regularized = n_randint(state, 3);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_add_ui(a1, a0, 1, prec0);

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_gamma_upper_asymp(w0, a0, z, regularized, prec0);
                break;
            case 1:
                acb_hypgeom_gamma_upper_1f1a(w0, a0, z, regularized, prec0);
                break;
            case 2:
                acb_hypgeom_gamma_upper_1f1b(w0, a0, z, regularized, prec0);
                break;
            default:
                acb_hypgeom_gamma_upper(w0, a0, z, regularized, prec0);
        }

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_gamma_upper_asymp(w1, a0, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_gamma_upper_1f1a(w1, a0, z, regularized, prec1);
                break;
            case 2:
                acb_hypgeom_gamma_upper_1f1b(w1, a0, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_gamma_upper(w1, a0, z, regularized, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 4))
        {
            case 0:
                acb_hypgeom_gamma_upper_asymp(w1, a1, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_gamma_upper_1f1a(w1, a1, z, regularized, prec1);
                break;
            case 2:
                acb_hypgeom_gamma_upper_1f1b(w1, a1, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_gamma_upper(w1, a1, z, regularized, prec1);
        }

        if (regularized == 2)
        {
            /* a R(a,z) + exp(-z) - z R(a+1,z) = 0 */
            acb_one(t);

            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_mul(t, t, u, prec0);

            acb_mul(b, w1, z, prec0);
            acb_addmul(t, a0, w0, prec0);
            acb_sub(t, t, b, prec0);
        }
        else if (regularized == 1)
        {
            /* Q(a,z) + exp(-z) z^a / Gamma(a+1) - Q(a+1,z) = 0 */
            /* http://dlmf.nist.gov/8.8.E6 */
            acb_pow(t, z, a0, prec0);
            acb_rgamma(u, a1, prec0);
            acb_mul(t, t, u, prec0);

            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_mul(t, t, u, prec0);

            acb_add(t, t, w0, prec0);
            acb_sub(t, t, w1, prec0);
        }
        else
        {
            /* a Gamma(a,z) + exp(-z) z^a - Gamma(a+1,z) = 0 */
            /* http://dlmf.nist.gov/8.8.E2 */
            acb_pow(t, z, a0, prec0);

            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_mul(t, t, u, prec0);

            acb_addmul(t, a0, w0, prec0);
            acb_sub(t, t, w1, prec0);
        }

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(t);
        acb_clear(u);
    }

    /* Accuracy regression tests. */
    {
        acb_t s, z;
        slong prec, issue, accuracy;
        acb_init(s);
        acb_init(z);

        issue = 166;
        prec = 165;
        accuracy = 100;
        acb_zero(s);
        acb_set_si(z, 110);
        _accuracy_regression_test(s, z, 2, prec, issue, accuracy);

        issue = 276;
        prec = 300;
        accuracy = 100;
        acb_set_ui(s, 357);
        acb_set_ui(z, 356);
        _accuracy_regression_test(s, z, 0, prec, issue, accuracy);
        arb_set_str(acb_realref(s), "356.123", prec);
        arb_set_str(acb_realref(z), "356.456", prec);
        _accuracy_regression_test(s, z, 0, prec, issue, accuracy);
        arb_set_str(acb_realref(s), "357.123", prec);
        arb_set_str(acb_realref(z), "356.456", prec);
        _accuracy_regression_test(s, z, 0, prec, issue, accuracy);
        arb_set_str(acb_realref(s), "357.456", prec);
        arb_set_str(acb_realref(z), "356.123", prec);
        _accuracy_regression_test(s, z, 0, prec, issue, accuracy);

        acb_clear(s);
        acb_clear(z);
    }

    /* Norm comparison tests (compare a^n to b^n + c^n). */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong prec;
        ulong n;
        arb_t a, b, c, u, v, w, rhs;
        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(u);
        arb_init(v);
        arb_init(w);
        arb_init(rhs);
        prec = n_randint(state, 1000) + 1;
        while (!arb_is_positive(a))
        {
            arb_randtest(a, state, n_randint(state, 1000)+1,
                                   n_randint(state, 100)+1);
        }
        while (!arb_is_positive(b))
        {
            arb_randtest(b, state, n_randint(state, 1000)+1,
                                   n_randint(state, 100)+1);
        }
        while (!arb_is_positive(c))
        {
            arb_randtest(c, state, n_randint(state, 1000)+1,
                                   n_randint(state, 100)+1);
        }
        if (n_randint(state, 20)) arb_zero(a);
        if (n_randint(state, 20)) arb_zero(b);
        if (n_randint(state, 20)) arb_zero(c);
        if (n_randint(state, 20)) arb_set(b, a);
        if (n_randint(state, 20)) arb_set(c, b);
        if (n_randint(state, 20)) arb_set(c, a);
        n = n_randint(state, 10);
        if (!n) n = WORD_MAX;
        if (n == WORD_MAX)
        {
            arb_set(u, a);
            arb_max(rhs, b, c, prec);
        }
        else
        {
            arb_pow_ui(u, a, n, prec);
            arb_pow_ui(v, b, n, prec);
            arb_pow_ui(w, c, n, prec);
            arb_add(rhs, v, w, prec);
        }
        if (arb_lt(u, rhs) || (arb_is_exact(u) && arb_equal(u, rhs)))
        {
            mag_t ma, mb, mc;
            mag_init(ma);
            mag_init(mb);
            mag_init(mc);
            arb_get_mag_lower(ma, a);
            arb_get_mag(mb, b);
            arb_get_mag(mc, c);
            if (_mag_gt_norm_ui(ma, mb, mc, n))
            {
                flint_printf("FAIL: _mag_gt_norm_ui\n\n");
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 30); flint_printf("\n\n");
                flint_printf("n = %ld\n\n", n);
                flint_abort();
            }
            mag_clear(ma);
            mag_clear(mb);
            mag_clear(mc);
        }
        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(u);
        arb_clear(v);
        arb_clear(w);
        arb_clear(rhs);
    }

    TEST_FUNCTION_END(state);
}

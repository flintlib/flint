/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_beta_lower, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, a1, b, b1, z, w, wa1, wb1, t, u;
        slong prec;
        int regularized;

        acb_init(a);
        acb_init(a1);
        acb_init(b);
        acb_init(b1);
        acb_init(z);
        acb_init(w);
        acb_init(wa1);
        acb_init(wb1);
        acb_init(t);
        acb_init(u);

        regularized = n_randint(state, 2);

        prec = 2 + n_randint(state, 100);

        acb_randtest_param(a, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_randtest_param(b, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));
        acb_randtest_param(z, state, 1 + n_randint(state, 100), 1 + n_randint(state, 10));

        acb_add_ui(a1, a, 1, prec);
        acb_add_ui(b1, b, 1, prec);

        acb_hypgeom_beta_lower(w, a, b, z, regularized, prec);
        acb_hypgeom_beta_lower(wa1, a1, b, z, regularized, prec);
        acb_hypgeom_beta_lower(wb1, a, b1, z, regularized, prec);

        if (regularized)
        {
            acb_mul(t, wa1, a, prec);
            acb_addmul(t, wb1, b, prec);
            acb_add(u, a, b, prec);
            acb_div(t, t, u, prec);
        }
        else
        {
            acb_add(t, wa1, wb1, prec);
        }

        if (!acb_overlaps(w, t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("wa1 = "); acb_printd(wa1, 30); flint_printf("\n\n");
            flint_printf("wb1 = "); acb_printd(wb1, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* test I(a,b;z) = 1-I(b,a,1-z) */

        if (!regularized)
        {
            acb_add(t, a, b, prec);
            acb_gamma(t, t, prec);
            acb_mul(w, w, t, prec);
            acb_rgamma(t, a, prec);
            acb_mul(w, w, t, prec);
            acb_rgamma(t, b, prec);
            acb_mul(w, w, t, prec);
        }

        acb_sub_ui(t, z, 1, prec);
        acb_neg(t, t);
        acb_hypgeom_beta_lower(t, b, a, t, 1, prec);
        acb_sub_ui(t, t, 1, prec);
        acb_neg(t, t);

        if (!acb_overlaps(w, t))
        {
            flint_printf("FAIL: symmetry\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(b1);
        acb_clear(z);
        acb_clear(w);
        acb_clear(wa1);
        acb_clear(wb1);
        acb_clear(t);
        acb_clear(u);
    }

    /* see issue 359 */
    {
        slong prec;
        acb_t z, a, b, res, res2;

        acb_init(z);
        acb_init(a);
        acb_init(b);
        acb_init(res);
        acb_init(res2);

        prec = 64;
        acb_set_d(a, 1e5);
        acb_set_d(b, 1e5);
        acb_set_ui(z, 4999);
        acb_div_ui(z, z, 10000, prec);
        acb_hypgeom_beta_lower(res, a, b, z, 1, prec);
        arb_set_str(acb_realref(res2), "0.4643650813520 +/- 5.17e-14", prec);

        if (!acb_overlaps(res, res2) || acb_rel_accuracy_bits(res) < acb_rel_accuracy_bits(res2) - 10)
        {
            flint_printf("FAIL: test case (1)\n\n");
            flint_printf("res = "); acb_printd(res, 100); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 100); flint_printf("\n\n");
            flint_abort();
        }

        prec = 128;
        acb_set_d(a, 1e10);
        acb_set_d(b, 1e10);
        acb_set_ui(z, 4999);
        acb_div_ui(z, z, 10000, prec);
        acb_hypgeom_beta_lower(res, a, b, z, 1, prec);
        arb_set_str(acb_realref(res2), "2.69791122252793228610950163e-176 +/- 2.47e-203", prec);

        if (!acb_overlaps(res, res2) || acb_rel_accuracy_bits(res) < acb_rel_accuracy_bits(res2) - 10)
        {
            flint_printf("FAIL: test case (2)\n\n");
            flint_printf("res = "); acb_printd(res, 100); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 100); flint_printf("\n\n");
            flint_abort();
        }

        prec = 128;
        acb_set_d(a, 1e15);
        acb_set_d(b, 1e15);
        acb_set_ui(z, 4999);
        acb_div_ui(z, z, 10000, prec);
        acb_hypgeom_beta_lower(res, a, b, z, 1, prec);
        arb_set_str(acb_realref(res2), "1.0612052723416047758478e-17371784 +/- 7.21e-17371807", prec);

        if (!acb_overlaps(res, res2) || acb_rel_accuracy_bits(res) < acb_rel_accuracy_bits(res2) - 10)
        {
            flint_printf("FAIL: test case (3\n\n");
            flint_printf("res = "); acb_printd(res, 100); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 100); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(z);
        acb_clear(a);
        acb_clear(b);
        acb_clear(res);
        acb_clear(res2);
    }

    TEST_FUNCTION_END(state);
}

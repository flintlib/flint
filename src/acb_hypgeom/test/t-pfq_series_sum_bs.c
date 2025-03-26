/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_pfq_series_sum_bs, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_struct a[4], b[4];
        acb_poly_t z, s1, s2, t1, t2;
        slong i, p, q, len1, len2, n, prec1, prec2;
        int regularized;

        p = n_randint(state, 3);
        q = n_randint(state, 3);

        len1 = n_randint(state, 8);
        len2 = n_randint(state, 8);

        prec1 = 2 + n_randint(state, 400);
        prec2 = 2 + n_randint(state, 400);

        n = n_randint(state, 50);

        regularized = n_randint(state, 2);

        acb_poly_init(z);
        acb_poly_init(s1);
        acb_poly_init(s2);
        acb_poly_init(t1);
        acb_poly_init(t2);
        for (i = 0; i < p; i++)
            acb_poly_init(a + i);
        for (i = 0; i < q; i++)
            acb_poly_init(b + i);

        acb_poly_randtest(z, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 10);

        for (i = 0; i < p; i++)
            acb_poly_randtest(a + i, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 3);
        for (i = 0; i < q; i++)
            acb_poly_randtest(b + i, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 3);

        acb_hypgeom_pfq_series_sum_forward(s1, t1, a, p, b, q, z, regularized, n, len1, prec1);
        acb_hypgeom_pfq_series_sum_bs(s2, t2, a, p, b, q, z, regularized, n, len2, prec2);

        acb_poly_truncate(s1, FLINT_MIN(len1, len2));
        acb_poly_truncate(s2, FLINT_MIN(len1, len2));
        acb_poly_truncate(t1, FLINT_MIN(len1, len2));
        acb_poly_truncate(t2, FLINT_MIN(len1, len2));

        if (!acb_poly_overlaps(s1, s2) || !acb_poly_overlaps(t1, t2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd\n", iter);
            flint_printf("n = %wd, len1 = %wd, len2 = %wd, prec1 = %wd, prec2 = %wd\n\n", n, len1, len2, prec1, prec2);
            flint_printf("p = %wd, q = %wd\n\n", p, q);
            flint_printf("z = "); acb_poly_printd(z, 15); flint_printf("\n\n");

            for (i = 0; i < p; i++)
            {
                flint_printf("a[%wd] = ", i); acb_poly_printd(a + i, 15); flint_printf("\n\n");
            }

            for (i = 0; i < q; i++)
            {
                flint_printf("b[%wd] = ", i); acb_poly_printd(b + i, 15); flint_printf("\n\n");
            }

            flint_printf("s1 = "); acb_poly_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); acb_poly_printd(s2, 15); flint_printf("\n\n");
            acb_poly_sub(s1, s1, s2, prec1);
            flint_printf("diff = "); acb_poly_printd(s1, 15); flint_printf("\n\n");

            flint_printf("t1 = "); acb_poly_printd(t1, 15); flint_printf("\n\n");
            flint_printf("t2 = "); acb_poly_printd(t2, 15); flint_printf("\n\n");
            acb_poly_sub(t1, t1, t2, prec1);
            flint_printf("diff = "); acb_poly_printd(t1, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_poly_clear(z);
        acb_poly_clear(s1);
        acb_poly_clear(s2);
        for (i = 0; i < p; i++)
            acb_poly_clear(a + i);
        for (i = 0; i < q; i++)
            acb_poly_clear(b + i);
        acb_poly_clear(t1);
        acb_poly_clear(t2);
    }


    /* issue #2268 */
    arb_t result;
    arb_init(result);

    slong p = 64;  // Precision

    // Define hypergeometric parameters: pFq({-10, 11, -10}, {1, -10}, 1)
    arb_struct *a, *b, z;
    slong p_len = 3, q_len = 2;  // 3F2

    a = _arb_vec_init(p_len);
    b = _arb_vec_init(q_len);
    arb_init(&z);

    // Set parameters
    arb_set_si(&a[0], -10);
    arb_set_si(&a[1], 11);
    arb_set_si(&a[2], -10);

    arb_set_si(&b[0], 1);
    arb_set_si(&b[1], -10);

    arb_set_si(&z, 1);  // z = 1

    // Compute the hypergeometric function
    arb_hypgeom_pfq(result, a, p_len, b, q_len, &z, 0, p);

    if (!arb_is_one(result))
    {
        flint_printf("FAIL: wrong result\n\n");
        flint_printf("HypergeometricPFQ[{-10, 11, -10}, {1, -10}, 1] = ");
        arb_print(result);
        flint_printf("\n");
        flint_abort();
    }

    // Clean up
    _arb_vec_clear(a, p_len);
    _arb_vec_clear(b, q_len);
    arb_clear(&z);
    arb_clear(result);

    TEST_FUNCTION_END(state);
}

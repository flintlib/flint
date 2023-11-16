/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_airy_zero, state)
{
    slong iter;

    /* interlacing test */
    {
        arb_t a, ap, b, bp, ap1;
        fmpz_t n;

        arb_init(a);
        arb_init(ap);
        arb_init(b);
        arb_init(bp);
        arb_init(ap1);
        fmpz_init(n);

        for (fmpz_one(n); fmpz_cmp_ui(n, 200) <= 0; fmpz_add_ui(n, n, 1))
        {
            arb_hypgeom_airy_zero(a, ap, b, bp, n, 53);
            fmpz_add_ui(n, n, 1);
            arb_hypgeom_airy_zero(NULL, ap1, NULL, NULL, n, 53);
            fmpz_sub_ui(n, n, 1);

            if (!arb_gt(ap, b) || !arb_gt(b, bp) || !arb_gt(bp, a) || !arb_gt(a, ap1))
            {
                flint_printf("FAIL: interlacing\n\n");
                flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
                flint_printf("a = "); arb_printn(a, 100, 0); flint_printf("\n\n");
                flint_printf("ap = "); arb_printn(ap, 100, 0); flint_printf("\n\n");
                flint_printf("b = "); arb_printn(b, 100, 0); flint_printf("\n\n");
                flint_printf("bp = "); arb_printn(bp, 100, 0); flint_printf("\n\n");
                flint_printf("ap1 = "); arb_printn(ap1, 100, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        for (fmpz_set_ui(n, 1000000); fmpz_cmp_ui(n, 1000000 + 200) <= 0; fmpz_add_ui(n, n, 1))
        {
            arb_hypgeom_airy_zero(a, ap, b, bp, n, 53);
            fmpz_add_ui(n, n, 1);
            arb_hypgeom_airy_zero(NULL, ap1, NULL, NULL, n, 53);
            fmpz_sub_ui(n, n, 1);

            if (!arb_gt(ap, b) || !arb_gt(b, bp) || !arb_gt(bp, a) || !arb_gt(a, ap1))
            {
                flint_printf("FAIL: interlacing\n\n");
                flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
                flint_printf("a = "); arb_printn(a, 100, 0); flint_printf("\n\n");
                flint_printf("ap = "); arb_printn(ap, 100, 0); flint_printf("\n\n");
                flint_printf("b = "); arb_printn(b, 100, 0); flint_printf("\n\n");
                flint_printf("bp = "); arb_printn(bp, 100, 0); flint_printf("\n\n");
                flint_printf("ap1 = "); arb_printn(ap1, 100, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(ap);
        arb_clear(b);
        arb_clear(bp);
        arb_clear(ap1);
        fmpz_clear(n);
    }

    /* self-consistency and accuracy test */
    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x1, x2, v1, v2;
        fmpz_t n;
        slong prec1, prec2;
        int which;

        arb_init(x1);
        arb_init(x2);
        arb_init(v1);
        arb_init(v2);
        fmpz_init(n);

        fmpz_randtest(n, state, 200);
        fmpz_abs(n, n);
        fmpz_add_ui(n, n, 1);
        prec1 = 2 + n_randtest(state) % 1000;
        prec2 = 2 + n_randtest(state) % 1000;
        which = n_randint(state, 4);

        if (which == 0)
        {
            arb_hypgeom_airy_zero(x1, NULL, NULL, NULL, n, prec1);
            arb_hypgeom_airy_zero(x2, NULL, NULL, NULL, n, prec2);
            arb_hypgeom_airy(v1, NULL, NULL, NULL, x1, prec1 + 30);
            arb_hypgeom_airy(v2, NULL, NULL, NULL, x2, prec2 + 30);
        }
        else if (which == 1)
        {
            arb_hypgeom_airy_zero(NULL, x1, NULL, NULL, n, prec1);
            arb_hypgeom_airy_zero(NULL, x2, NULL, NULL, n, prec2);
            arb_hypgeom_airy(NULL, v1, NULL, NULL, x1, prec1 + 30);
            arb_hypgeom_airy(NULL, v2, NULL, NULL, x2, prec2 + 30);
        }
        else if (which == 2)
        {
            arb_hypgeom_airy_zero(NULL, NULL, x1, NULL, n, prec1);
            arb_hypgeom_airy_zero(NULL, NULL, x2, NULL, n, prec2);
            arb_hypgeom_airy(NULL, NULL, v1, NULL, x1, prec1 + 30);
            arb_hypgeom_airy(NULL, NULL, v2, NULL, x2, prec2 + 30);
        }
        else
        {
            arb_hypgeom_airy_zero(NULL, NULL, NULL, x1, n, prec1);
            arb_hypgeom_airy_zero(NULL, NULL, NULL, x2, n, prec2);
            arb_hypgeom_airy(NULL, NULL, NULL, v1, x1, prec1 + 30);
            arb_hypgeom_airy(NULL, NULL, NULL, v2, x2, prec2 + 30);
        }

        if (!arb_overlaps(x1, x2) || !arb_contains_zero(v1) || !arb_contains_zero(v2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("which = %d, n = ", which); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); arb_printn(x1, 100, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 100, 0); flint_printf("\n\n");
            flint_printf("v1 = "); arb_printn(v1, 100, 0); flint_printf("\n\n");
            flint_printf("v2 = "); arb_printn(v2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(x1) < prec1 - 3 || arb_rel_accuracy_bits(x2) < prec2 - 3)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("which = %d, n = ", which); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("acc(x1) = %wd, acc(x2) = %wd\n\n", arb_rel_accuracy_bits(x1), arb_rel_accuracy_bits(x2));
            flint_printf("x1 = "); arb_printn(x1, 100, 0); flint_printf("\n\n");
            flint_printf("x2 = "); arb_printn(x2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x1);
        arb_clear(x2);
        arb_clear(v1);
        arb_clear(v2);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}

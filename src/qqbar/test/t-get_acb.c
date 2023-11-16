/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_get_acb, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x;
        acb_t z, w;
        slong prec1, prec2, acc1, acc2;

        qqbar_init(x);
        acb_init(z);
        acb_init(w);

        qqbar_randtest(x, state, 3, 400);
        prec1 = 2 + n_randint(state, 800);
        prec2 = prec1 + 2 + n_randint(state, 400);

        qqbar_get_acb(z, x, prec1);
        qqbar_get_acb(w, x, prec2);

        acc1 = arb_rel_accuracy_bits(acb_realref(z));
        acc2 = arb_rel_accuracy_bits(acb_imagref(z));

        if (acc1 < prec1 - 2 || acc2 < prec1 - 2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("z = "); acb_printn(z, 1000, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec1, acc1, acc2);
            flint_abort();
        }

        acc1 = arb_rel_accuracy_bits(acb_realref(w));
        acc2 = arb_rel_accuracy_bits(acb_imagref(w));

        if (acc1 < prec2 - 2 || acc2 < prec2 - 2)
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("w = "); acb_printn(w, 1000, ARB_STR_CONDENSE * 30); flint_printf("\n\n");
            flint_printf("prec2 = %wd, acc1 = %wd, acc2 = %wd\n\n", prec2, acc1, acc2);
            flint_abort();
        }

        if (!acb_overlaps(z, w))
        {
            flint_printf("FAIL!\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("z re = "); arb_printn(acb_realref(z), 1000, 0); flint_printf("\n\n");
            flint_printf("z im = "); arb_printn(acb_imagref(z), 1000, 0); flint_printf("\n\n");
            flint_printf("w re = "); arb_printn(acb_realref(w), 1000, 0); flint_printf("\n\n");
            flint_printf("w im = "); arb_printn(acb_imagref(w), 1000, 0); flint_printf("\n\n");
            flint_printf("z re = "); arb_print(acb_realref(z)); flint_printf("\n\n");
            flint_printf("z im = "); arb_print(acb_imagref(z)); flint_printf("\n\n");
            flint_printf("w re = "); arb_print(acb_realref(w)); flint_printf("\n\n");
            flint_printf("w im = "); arb_print(acb_imagref(w)); flint_printf("\n\n");
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_abort();
        }

        qqbar_clear(x);
        acb_clear(z);
        acb_clear(w);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_lambertw, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x1, x2, t, w1, w2;
        slong prec1, prec2, ebits;
        int flags;
        fmpz_t k;

        acb_init(x1);
        acb_init(x2);
        acb_init(t);
        acb_init(w1);
        acb_init(w2);
        fmpz_init(k);

        if (n_randint(state, 4) == 0)
        {
            prec1 = 2 + n_randint(state, 3000);
            prec2 = 2 + n_randint(state, 3000);
            ebits = 1 + n_randint(state, 1000);
        }
        else
        {
            prec1 = 2 + n_randint(state, 300);
            prec2 = 2 + n_randint(state, 300);
            ebits = 1 + n_randint(state, 50);
        }

        acb_randtest(x1, state, 1 + n_randint(state, 1000), ebits);
        acb_randtest(x2, state, 1 + n_randint(state, 1000), ebits);
        acb_randtest(t, state, 1 + n_randint(state, 1000), ebits);
        acb_randtest(w1, state, 1 + n_randint(state, 1000), ebits);
        acb_randtest(w2, state, 1 + n_randint(state, 1000), ebits);
        fmpz_randtest(k, state, ebits);

        flags = 0;

        switch (n_randint(state, 10))
        {
            case 0:
                flags = ACB_LAMBERTW_LEFT;
                break;
            case 1:
                flags = ACB_LAMBERTW_MIDDLE;
                fmpz_set_si(k, -1);
                break;
            default:
                break;
        }

        if (n_randint(state, 4) == 0)
        {
            arb_const_e(acb_realref(t), 2 * prec1);
            arb_inv(acb_realref(t), acb_realref(t), 2 * prec1);
            arb_sub(acb_realref(x1), acb_realref(x1), acb_realref(t), 2 * prec1);
        }

        if (n_randint(state, 2))
        {
            acb_set(x2, x1);
        }
        else
        {
            acb_add(x2, x1, t, 2 * prec1);
            acb_sub(x2, x2, t, 2 * prec1);
        }

        if (n_randint(state, 4) == 0 && flags == 0)
            acb_lambertw_asymp(w1, x1, k,
                1 + n_randint(state, 10), 1 + n_randint(state, 10), prec1);
        else
            acb_lambertw(w1, x1, k, flags, prec1);

        if (n_randint(state, 2))
        {
            acb_set(w2, x2);
            acb_lambertw(w2, w2, k, flags, prec1);
        }
        else
        {
            acb_lambertw(w2, x2, k, flags, prec1);
        }

        if (!acb_overlaps(w1, w2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter %wd, flags = %d, branch = ", flags, iter); fmpz_print(k);
            flint_printf(" prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); acb_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); acb_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
            flint_abort();
        }

        acb_exp(t, w1, prec1);
        acb_mul(t, t, w1, prec1);

        if (!acb_contains(t, x1))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("iter %wd, flags = %d, branch = ", flags, iter); fmpz_print(k);
            flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
            flint_printf("x1 = "); acb_printd(x1, 50); flint_printf("\n\n");
            flint_printf("x2 = "); acb_printd(x2, 50); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 50); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
            flint_printf("t  = "); acb_printd(t, 50); flint_printf("\n\n");
            flint_abort();
        }

        if (!arb_contains_zero(acb_imagref(x1)))
        {
            acb_conj(x2, x1);
            fmpz_neg(k, k);
            if (flags == 2)
                fmpz_sub_ui(k, k, 1);

            acb_lambertw(w2, x2, k, flags, prec2);

            if (flags == 2)
                fmpz_add_ui(k, k, 1);
            fmpz_neg(k, k);
            acb_conj(w2, w2);

            if (!acb_overlaps(w1, w2))
            {
                flint_printf("FAIL: conjugation\n\n");
                flint_printf("iter %wd, flags = %d, branch = ", flags, iter); fmpz_print(k);
                flint_printf("prec1 = %wd, prec2 = %wd\n\n", prec1, prec2);
                flint_printf("x1 = "); acb_printd(x1, 50); flint_printf("\n\n");
                flint_printf("x2 = "); acb_printd(x2, 50); flint_printf("\n\n");
                flint_printf("w1 = "); acb_printd(w1, 50); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(x1);
        acb_clear(x2);
        acb_clear(t);
        acb_clear(w1);
        acb_clear(w2);
        fmpz_clear(k);
    }

    TEST_FUNCTION_END(state);
}

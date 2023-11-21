/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "fexpr.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_get_fexpr_formula, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        fexpr_t e;
        ulong flags;

        qqbar_init(x);
        qqbar_init(y);
        fexpr_init(e);

        qqbar_randtest(x, state, 5, 10);

        if (n_randint(state, 2))
            flags = n_randlimb(state);
        else
            flags = QQBAR_FORMULA_ALL;

        if (qqbar_get_fexpr_formula(e, x, flags))
        {
            qqbar_set_fexpr(y, e);

            if (!qqbar_equal(x, y))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
                flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
                flint_abort();
            }
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fexpr_clear(e);
    }

    /* test trigonometrics and exponentials */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        qqbar_t x, y;
        fexpr_t e;
        slong p, q;
        fmpq_t a;

        qqbar_init(x);
        qqbar_init(y);
        fexpr_init(e);
        fmpq_init(a);

        p = -10 + n_randint(state, 100);

        if (n_randint(state, 10))
            q = 1 + n_randint(state, 6);
        else
            q = 1 + n_randint(state, 10);

        switch (n_randint(state, 7))
        {
            case 0: qqbar_exp_pi_i(x, p, q); break;
            case 1: qqbar_sin_pi(x, p, q); break;
            case 2: qqbar_cos_pi(x, p, q); break;
            case 3: qqbar_tan_pi(x, p, q); break;
            case 4: qqbar_cot_pi(x, p, q); break;
            case 5: qqbar_sec_pi(x, p, q); break;
            case 6: qqbar_csc_pi(x, p, q); break;
        }

        fmpq_randtest(a, state, 3);
        qqbar_add_fmpq(x, x, a);
        fmpq_randtest(a, state, 3);
        qqbar_mul_fmpq(x, x, a);

        if (qqbar_degree(x) <= 2)
        {
            p = -10 + n_randint(state, 100);
            q = 1 + n_randint(state, 6);

            switch (n_randint(state, 7))
            {
                case 0: qqbar_exp_pi_i(y, p, q); break;
                case 1: qqbar_sin_pi(y, p, q); break;
                case 2: qqbar_cos_pi(y, p, q); break;
                case 3: qqbar_tan_pi(y, p, q); break;
                case 4: qqbar_cot_pi(y, p, q); break;
                case 5: qqbar_sec_pi(y, p, q); break;
                case 6: qqbar_csc_pi(y, p, q); break;
            }

            qqbar_add(x, x, y);
            qqbar_zero(y);
        }

        if (qqbar_get_fexpr_formula(e, x, QQBAR_FORMULA_ALL))
        {
            qqbar_set_fexpr(y, e);

            if (!qqbar_equal(x, y))
            {
                flint_printf("FAIL!\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
                flint_printf("e = "); fexpr_print(e); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            flint_printf("Unexpected:\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
        }

        qqbar_clear(x);
        qqbar_clear(y);
        fexpr_clear(e);
        fmpq_clear(a);
    }

    TEST_FUNCTION_END(state);
}

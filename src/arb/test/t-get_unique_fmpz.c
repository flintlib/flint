/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_get_unique_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, z;
        fmpz_t y, a, b, exp;
        int unique, unique2;

        arb_init(x);
        arb_init(z);
        fmpz_init(y);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(exp);

        arb_randtest(x, state, 2000, 10);

        /* generate tiny and huge radii */
        if (iter % 2 == 0)
        {
            mag_randtest_special(arb_radref(x), state, 100);

            unique = arb_get_unique_fmpz(y, x);

            arf_get_fmpz(a, arb_midref(x), ARF_RND_FLOOR);
            fmpz_add_ui(b, a, 1);

            if (unique)
            {
                if (arb_contains_fmpz(x, a) == arb_contains_fmpz(x, b))
                {
                    flint_printf("FAIL (1):\n\n");
                    flint_printf("x = "); arb_printd(x, 100); flint_printf("\n\n");
                    flint_printf("unique = %d\n\n", unique);
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
                    flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
                    flint_abort();
                }
            }
            else
            {
                if (arb_contains_fmpz(x, a) != arb_contains_fmpz(x, b))
                {
                    flint_printf("FAIL (2):\n\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("unique = %d\n\n", unique);
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
                    flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }
        else
        {
            unique = arb_get_unique_fmpz(y, x);

            arb_get_interval_fmpz_2exp(a, b, exp, x);
            if (fmpz_sgn(exp) >= 0)
            {
                fmpz_mul_2exp(a, a, fmpz_get_si(exp));
                fmpz_mul_2exp(b, b, fmpz_get_si(exp));
            }
            else
            {
                fmpz_cdiv_q_2exp(a, a, -fmpz_get_si(exp));
                fmpz_fdiv_q_2exp(b, b, -fmpz_get_si(exp));
            }
            unique2 = fmpz_equal(a, b);

            if ((unique != unique2) || (unique && !fmpz_equal(y, a)))
            {
                flint_printf("FAIL:\n\n");
                flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                flint_printf("unique = %d, unique2 = %d\n\n", unique, unique2);
                flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
                flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
                flint_printf(" exp = "); fmpz_print(exp); flint_printf("\n\n");
                flint_abort();
            }

            if (unique)
            {
                arb_set_fmpz(z, y);
                arb_set_round(z, z, 2 + n_randint(state, 1000));

                if (!arb_overlaps(x, z))
                {
                    flint_printf("FAIL (overlap):\n\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arb_print(z); flint_printf("\n\n");
                    flint_abort();
                }

                fmpz_add_ui(b, y, 1);
                if (arb_contains_fmpz(x, b))
                {
                    flint_printf("FAIL (contains a + 1):\n\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_abort();
                }

                fmpz_sub_ui(b, y, 1);
                if (arb_contains_fmpz(x, b))
                {
                    flint_printf("FAIL (contains a - 1):\n\n");
                    flint_printf("x = "); arb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        arb_clear(x);
        arb_clear(z);
        fmpz_clear(y);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(exp);
    }

    TEST_FUNCTION_END(state);
}

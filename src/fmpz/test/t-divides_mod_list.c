/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_divides_mod_list, state)
{
    slong i, j;
    int result1, result2;

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fmpz_t xstart, xstride, xlength;
        fmpz_t a, b, n;
        fmpz_t q, r, x;

        fmpz_init(xstart);
        fmpz_init(xstride);
        fmpz_init(xlength);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(n);
        fmpz_init(q);
        fmpz_init(r);
        fmpz_init(x);

        fmpz_randtest(xstart, state, 100);
        fmpz_randtest(xstride, state, 100);
        fmpz_randtest(xlength, state, 100);

        fmpz_randtest(a, state, 100);
        fmpz_randtest(b, state, 100);
        fmpz_randtest_unsigned(n, state, 100);
        fmpz_add_ui(n, n, 1);

        fmpz_divides_mod_list(xstart, xstride, xlength, a, b, n);

        if (fmpz_sgn(xlength) > 0)
        {
            if (fmpz_sgn(xstride) <= 0)
            {
                flint_printf("FAIL:\ncheck xstride > 0\n");
                fflush(stdout);
                flint_abort();
            }

            if (fmpz_cmp(xstart, xstride) >= 0)
            {
                flint_printf("FAIL:\ncheck xstart < xstride\n");
                fflush(stdout);
                flint_abort();
            }

            for (j = 0; j < 200; j++)
            {
                fmpz_randm(x, state, n);
                fmpz_mul(r, x, b);
                fmpz_sub(r, r, a);
                result1 = fmpz_divisible(r, n);

                fmpz_sub(r, x, xstart);
                fmpz_fdiv_qr(q, r, r, xstride);
                result2 = fmpz_is_zero(r) &&
                          fmpz_sgn(q) >= 0 &&
                          fmpz_cmp(q, xlength) < 0;

                if (result1 != result2)
                {
                    flint_printf("FAIL:\ncheck ");
                    fmpz_print(a);
                    flint_printf("/");
                    fmpz_print(b);
                    flint_printf(" = ");
                    fmpz_print(x);
                    flint_printf(" mod ");
                    fmpz_print(n);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
        else
        {
            for (j = 0; j < 200; j++)
            {
                fmpz_randm(x, state, n);
                fmpz_mul(r, x, b);
                fmpz_sub(r, r, a);
                result1 = fmpz_divisible(r, n);

                if (result1)
                {
                    flint_printf("FAIL:\ncheck ");
                    fmpz_print(a);
                    flint_printf("/");
                    fmpz_print(b);
                    flint_printf(" = ");
                    fmpz_print(x);
                    flint_printf(" mod ");
                    fmpz_print(n);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_set(q, a);
        fmpz_set(r, b);
        fmpz_set(x, n);
        fmpz_divides_mod_list(q, r, x, q, r, x);
        if (!fmpz_equal(q, xstart) ||
            !fmpz_equal(r, xstride) ||
            !fmpz_equal(x, xlength))
        {
            flint_printf("FAIL:\ncheck aliasing 1\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_set(x, a);
        fmpz_set(r, b);
        fmpz_set(q, n);
        fmpz_divides_mod_list(q, r, x, x, r, q);
        if (!fmpz_equal(q, xstart) ||
            !fmpz_equal(r, xstride) ||
            !fmpz_equal(x, xlength))
        {
            flint_printf("FAIL:\ncheck aliasing 2\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_set(r, a);
        fmpz_set(x, b);
        fmpz_set(q, n);
        fmpz_divides_mod_list(q, r, x, r, x, q);
        if (!fmpz_equal(q, xstart) ||
            !fmpz_equal(r, xstride) ||
            !fmpz_equal(x, xlength))
        {
            flint_printf("FAIL:\ncheck aliasing 3\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(xstart);
        fmpz_clear(xstride);
        fmpz_clear(xlength);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(n);
        fmpz_clear(q);
        fmpz_clear(r);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}

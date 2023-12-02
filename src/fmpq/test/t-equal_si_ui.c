/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_equal_si_ui, state)
{
    int i;
    int res;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_t x;
        slong y_si = 0; /* Silence warnings about maybe uninitialized */
        ulong y_ui = 0;
        int type;

        fmpq_init(x);

        type = n_randint(state, 4);

        if (type == 0)
        {
            /* Let x be proper rational, y is slong */
            do { fmpq_randtest(x, state, 200); } while (x->den == WORD(1));
            y_si = z_randtest(state);

            /* Then x != y */
            res = (fmpq_equal_si(x, y_si) == 0);
        }
        else if (type == 1)
        {
            /* Let x be integer, y is slong */
            fmpz_one(fmpq_denref(x));
            y_si = z_randtest(state);

            if (n_randint(state, 2))
            {
                /* Random numerator */
                fmpz_randtest(fmpq_numref(x), state, 200);
                res = (fmpq_equal_si(x, y_si) == fmpz_equal_si(fmpq_numref(x), y_si));
            }
            else
            {
                /* Numerator of x equals to y */
                fmpz_set_si(fmpq_numref(x), y_si);
                res = (fmpq_equal_si(x, y_si) == 1);
            }
        }
        else if (type == 2)
        {
            /* Let x be proper rational, y is ulong */
            do { fmpq_randtest(x, state, 200); } while (x->den == WORD(1));
            y_ui = n_randtest(state);

            /* Then x != y */
            res = (fmpq_equal_ui(x, y_ui) == 0);
        }
        else
        {
            /* Let x be integer, y is ulong */
            fmpz_one(fmpq_denref(x));
            y_ui = n_randtest(state);

            if (n_randint(state, 2))
            {
                /* Random numerator */
                fmpz_randtest(fmpq_numref(x), state, 200);
                res = (fmpq_equal_ui(x, y_ui) == fmpz_equal_ui(fmpq_numref(x), y_ui));
            }
            else
            {
                /* Numerator of x equals to y */
                fmpz_set_ui(fmpq_numref(x), y_ui);
                res = (fmpq_equal_ui(x, y_ui) == 1);
            }
        }

        if (!res)
        {
            flint_printf("FAIL\n");
            flint_printf("x = ");
            fmpq_print(x);
            flint_printf("\ntype = %d", type);
            if (type <= 1)
                flint_printf("y_si = %wd\n", y_si);
            else
                flint_printf("y_ui = %wu\n", y_ui);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
    }

    TEST_FUNCTION_END(state);
}

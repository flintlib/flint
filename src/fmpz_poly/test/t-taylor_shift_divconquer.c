/*
    Copyright (C) 2012, 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_taylor_shift_divconquer, state)
{
    int i;

    /* Compare with composition */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, h1, h2;
        fmpz_t c;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h1);
        fmpz_poly_init(h2);

        fmpz_init(c);

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 1000),
                                     1 + n_randint(state, 400));

        fmpz_randtest(c, state, n_randint(state, 10));

        fmpz_poly_set_coeff_ui(g, 1, 1);
        fmpz_poly_set_coeff_fmpz(g, 0, c);

        flint_set_num_threads(1 + n_randint(state, 3));

        if (n_randint(state, 2))
        {
            fmpz_poly_taylor_shift_divconquer(h1, f, c);
        }
        else
        {
            fmpz_poly_set(h1, f);
            fmpz_poly_taylor_shift_divconquer(h1, h1, c);
        }

        flint_set_num_threads(1);

        fmpz_poly_compose_divconquer(h2, f, g);

        if (!fmpz_poly_equal(h1, h2))
        {
            flint_printf("FAIL\n");
            fmpz_poly_print(f); flint_printf("\n");
            fmpz_poly_print(g); flint_printf("\n");
            fmpz_poly_print(h1); flint_printf("\n");
            fmpz_poly_print(h2); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h1);
        fmpz_poly_clear(h2);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}

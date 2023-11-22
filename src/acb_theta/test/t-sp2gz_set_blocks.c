/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_sp2gz_set_blocks, state)
{
    slong iter;

    /* Test: set_abcd is inverse of get_abcd */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        fmpz_mat_t a, b, c, d;
        fmpz_mat_t m, n;
        slong bits = n_randint(state, 10);

        fmpz_mat_init(m, 2 * g, 2 * g);
        fmpz_mat_init(n, 2 * g, 2 * g);
        sp2gz_randtest(m, state, bits);

        fmpz_mat_window_init(a, m, 0, 0, g, g);
        fmpz_mat_window_init(b, m, 0, g, g, 2 * g);
        fmpz_mat_window_init(c, m, g, 0, 2 * g, g);
        fmpz_mat_window_init(d, m, g, g, 2 * g, 2 * g);

        sp2gz_set_blocks(n, a, b, c, d);

        if (!fmpz_mat_equal(m, n))
        {
            flint_printf("FAIL\n\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\n\n");
            fmpz_mat_print_pretty(n);
            flint_abort();
            flint_printf("\n\n");
        }

        fmpz_mat_window_clear(a);
        fmpz_mat_window_clear(b);
        fmpz_mat_window_clear(c);
        fmpz_mat_window_clear(d);
        fmpz_mat_clear(m);
        fmpz_mat_clear(n);
    }

    TEST_FUNCTION_END(state);
}

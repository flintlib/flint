/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_g2_character, state)
{
    slong iter;

    /* Test: agrees with kappa2 */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t mat;
        slong bits = n_randint(state, 10);
        ulong * chars;
        slong * es;
        slong ab, eps, test;

        fmpz_mat_init(mat, 4, 4);
        sp2gz_randtest(mat, state, bits);
        chars = flint_malloc(16 * sizeof(ulong));
        es = flint_malloc(16 * sizeof(slong));

        eps = acb_theta_g2_character(mat);

        test = 10 * acb_siegel_kappa2(mat); /* 10 theta constants */
        acb_theta_char_table(chars, es, mat, 0, 1);
        for (ab = 0; ab < 16; ab++)
        {
            if (acb_theta_char_is_even(ab, 2))
            {
                test += es[ab];
            }
        }
        if (test % 4 != 0)
        {
            flint_printf("FAIL (%wd mod 4)\n", test % 4);
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            flint_abort();
        }
        test = (test / 4) % 2;

        if (eps != test)
        {
            flint_printf("FAIL (%wd != %wd)\n", eps, test);
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(mat);
        flint_free(chars);
        flint_free(es);
    }

    TEST_FUNCTION_END(state);
}

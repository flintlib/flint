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

TEST_FUNCTION_START(acb_theta_char_table, state)
{
    slong iter;

    /* Test: on trigonal symplectic matrices, a remains the same */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n2 = 1 << (2 * g);
        slong bits = 8;
        fmpz_mat_t mat;
        ulong * chars;
        slong * es;
        ulong ab;
        slong j, k;

        fmpz_mat_init(mat, 2 * g, 2 * g);
        chars = flint_malloc(n2 * sizeof(ulong));
        es = flint_malloc(n2 * sizeof(slong));

        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                fmpz_randtest(fmpz_mat_entry(mat, j, k), state, bits);
                fmpz_set(fmpz_mat_entry(mat, k, j), fmpz_mat_entry(mat, j, k));
            }
        }
        sp2gz_trig(mat, mat);

        acb_theta_char_table(chars, es, mat, -1);
        for (ab = 0; ab < n2; ab++)
        {
            if ((chars[ab] >> g) != (ab >> g))
            {
                flint_printf("FAIL\n");
                flint_printf("ab = %wd, test = %wd, matrix:\n", ab, chars[ab]);
                fmpz_mat_print_pretty(mat);
                flint_abort();
            }
        }

        fmpz_mat_clear(mat);
        flint_free(chars);
        flint_free(es);
    }

    TEST_FUNCTION_END(state);
}

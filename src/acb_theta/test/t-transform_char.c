/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("transform_char....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: on trigonal symplectic matrices, a remains the same */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong bits = 8;
        fmpz_mat_t mat;
        slong e;
        ulong ab = n_randint(state, 1 << (2 * g));
        ulong test;
        slong j, k;

        fmpz_mat_init(mat, 2 * g, 2 * g);

        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                fmpz_randtest(fmpz_mat_entry(mat, j, k), state, bits);
                fmpz_set(fmpz_mat_entry(mat, k, j), fmpz_mat_entry(mat, j, k));
            }
        }
        sp2gz_trig(mat, mat);

        test = acb_theta_transform_char(&e, mat, ab);

        if ((test >> g) != (ab >> g))
        {
            flint_printf("FAIL\n");
            flint_printf("ab = %wd, test = %wd, matrix:\n", ab, test);
            fmpz_mat_print_pretty(mat);
            flint_abort();
        }

        fmpz_mat_clear(mat);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

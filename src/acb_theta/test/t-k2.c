/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("k2....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: value on [u, 0; 0, u^-t] is det(u) */
    for (iter = 0; iter < 50 * arb_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        fmpz_mat_t U, mat;
        fmpz_t det;
        slong k2;
        slong bits = 1 + n_randint(state, 10);

        fmpz_mat_init(U, g, g);
        fmpz_mat_init(mat, 2 * g, 2 * g);
        fmpz_init(det);

        fmpz_mat_one(U);
        if (iter % 2 == 0)
            fmpz_set_si(fmpz_mat_entry(U, 0, 0), -1);

        fmpz_mat_randops(U, state, 2 * bits);
        fmpz_mat_diag_sp(mat, U);
        fmpz_mat_det(det, U);
        k2 = acb_theta_k2(mat);

        /* det is 1 or -1; k2 is 0 or 2 */
        if (k2 != 1 - fmpz_get_si(det))
        {
            flint_printf("FAIL\n");
            fmpz_mat_print_pretty(mat);
            flint_printf("\n");
            flint_printf("k2: %wd\n", k2);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(U);
        fmpz_mat_clear(mat);
        fmpz_clear(det);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

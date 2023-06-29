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

    flint_printf("siegel_reduce_imag....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: is almost Minkowski reduction for g=2 */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 2;
        slong prec = 100 + n_randint(state, 500);
        slong mag_bits = 1 + n_randint(state, 5);
        arb_mat_t M;
        arb_mat_t R;
        arb_mat_t T;
        fmpz_mat_t U;
        arb_t test;
        int res;

        arb_mat_init(M, g, g);
        arb_mat_init(R, g, g);
        arb_mat_init(T, g, g);
        fmpz_mat_init(U, g, g);
        arb_init(test);

        arb_mat_randtest_spd(M, state, prec, mag_bits);
        arb_mat_spd_lll_reduce(U, M, prec);

        arb_mat_set_fmpz_mat(T, U);
        arb_mat_mul(R, T, M, prec);
        arb_mat_transpose(T, T);
        arb_mat_mul(R, R, T, prec);

        arb_set_d(test, 1.98);
        arb_mul(test, test, arb_mat_entry(R, 0, 1), prec);
        arb_abs(test, test);
        res = arb_gt(arb_mat_entry(R, 0, 0), test);

        arb_set_d(test, 1.02);
        arb_mul(test, test, arb_mat_entry(R, 1, 1), prec);
        res = res && arb_gt(test, arb_mat_entry(R, 0, 0));

        if (!res)
        {
            flint_printf("FAIL (reduction)\n");
            arb_mat_printd(M, 10);
            arb_mat_printd(R, 10);
            fmpz_mat_print_pretty(U);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        arb_mat_clear(M);
        arb_mat_clear(R);
        arb_mat_clear(T);
        fmpz_mat_clear(U);
        arb_clear(test);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

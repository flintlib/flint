/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("dct....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        arb_mat_t A, Ainv, AT;
        slong n, prec;

        n = n_randint(state, 10);
        prec = 53 + n_randint(state, 30);

        arb_mat_init(A, n, n);
        arb_mat_init(Ainv, n, n);
        arb_mat_init(AT, n, n);

        arb_mat_randtest(A, state, 100, 10);
        arb_mat_dct(A, 0, prec);

        if (!arb_mat_inv(Ainv, A, prec))
        {
            flint_printf("FAIL: small DCT matrix (n = %wd) not invertible\n", n);
            flint_abort();
        }

        arb_mat_transpose(AT, A);

        if (!arb_mat_overlaps(AT, Ainv))
        {
            flint_printf("FAIL: overlap (n = %wd)\n", n);
            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("Ainv = \n"); arb_mat_printd(Ainv, 15); flint_printf("\n\n");
            flint_printf("AT = \n"); arb_mat_printd(AT, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(Ainv);
        arb_mat_clear(AT);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

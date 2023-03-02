/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("dft....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        acb_mat_t A, Ainv, AT;
        slong n, prec;

        n = n_randint(state, 10);
        prec = 53 + n_randint(state, 30);

        acb_mat_init(A, n, n);
        acb_mat_init(Ainv, n, n);
        acb_mat_init(AT, n, n);

        acb_mat_randtest(A, state, 100, 10);
        acb_mat_dft(A, 0, prec);

        if (!acb_mat_inv(Ainv, A, prec))
        {
            flint_printf("FAIL: small DFT matrix (n = %wd) not invertible\n", n);
            flint_abort();
        }

        acb_mat_conjugate_transpose(AT, A);

        if (!acb_mat_overlaps(AT, Ainv))
        {
            flint_printf("FAIL: overlap (n = %wd)\n", n);
            flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
            flint_printf("Ainv = \n"); acb_mat_printd(Ainv, 15); flint_printf("\n\n");
            flint_printf("AT = \n"); acb_mat_printd(AT, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_mat_clear(A);
        acb_mat_clear(Ainv);
        acb_mat_clear(AT);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

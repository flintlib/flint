/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mat.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);


    flint_printf("lll_original....");
    fflush(stdout);

    /* check output basis is LLL reduced (randajtai used) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int result;
        fmpz_mat_t A;
        fmpq_t delta, eta;

        slong m;

        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpq_init(delta);
        fmpq_init(eta);

        fmpq_set_si(delta, 3, 4);
        fmpq_set_si(eta, 1, 2);

        fmpz_mat_randajtai(A, state, 0.5);

        fmpz_mat_lll_original(A, delta, eta);

        result = fmpz_mat_is_reduced(A, 0.75, 0.5);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            fmpz_mat_print_pretty(A);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpq_clear(delta);
        fmpq_clear(eta);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

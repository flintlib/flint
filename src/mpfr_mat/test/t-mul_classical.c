/*
    Copyright (C) 2010 Fredrik Johansson
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
#include "mpfr_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("mul_classical....");
    fflush(stdout);

    /* Check aliasing C and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mpfr_mat_t A, B, C;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        mpfr_mat_init(A, m, n, 200);
        mpfr_mat_init(B, n, n, 200);
        mpfr_mat_init(C, m, n, 200);

        mpfr_mat_randtest(A, state);
        mpfr_mat_randtest(B, state);
        mpfr_mat_randtest(C, state);

        mpfr_mat_mul_classical(C, A, B, MPFR_RNDN);
        mpfr_mat_mul_classical(A, A, B, MPFR_RNDN);

        if (!mpfr_mat_equal(C, A))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        mpfr_mat_clear(A);
        mpfr_mat_clear(B);
        mpfr_mat_clear(C);
    }

    /* Check aliasing C and B */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        mpfr_mat_t A, B, C;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        mpfr_mat_init(A, m, m, 200);
        mpfr_mat_init(B, m, n, 200);
        mpfr_mat_init(C, m, n, 200);

        mpfr_mat_randtest(A, state);
        mpfr_mat_randtest(B, state);
        mpfr_mat_randtest(C, state);

        mpfr_mat_mul_classical(C, A, B, MPFR_RNDN);
        mpfr_mat_mul_classical(B, A, B, MPFR_RNDN);

        if (!mpfr_mat_equal(C, B))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        mpfr_mat_clear(A);
        mpfr_mat_clear(B);
        mpfr_mat_clear(C);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

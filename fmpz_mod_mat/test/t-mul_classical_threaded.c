/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_mat.h"
#include "ulong_extras.h"
#include "thread_support.h"

int
main(void)
{
#if HAVE_PTHREAD && (HAVE_TLS || FLINT_REENTRANT)
    slong i, max_threads = 5;
#endif
    FLINT_TEST_INIT(state);

    flint_printf("mul_classical_threaded....");
    fflush(stdout);

#if HAVE_PTHREAD && (HAVE_TLS || FLINT_REENTRANT)
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, C, D;
        fmpz_t mod;

        slong m, k, n;

        fmpz_init(mod);

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        fmpz_randtest_not_zero(mod, state, 100);

        fmpz_mod_mat_init(A, m, n, mod);
        fmpz_mod_mat_init(B, n, k, mod);
        fmpz_mod_mat_init(C, m, k, mod);
        fmpz_mod_mat_init(D, m, k, mod);

        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);
        fmpz_mod_mat_randtest(C, state);  /* make sure noise in the output is ok */

        fmpz_mod_mat_mul_classical_threaded(C, A, B);
        fmpz_mod_mat_mul(D, A, B);

        if (!fmpz_mod_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            fmpz_mod_mat_print_pretty(A);
            fmpz_mod_mat_print_pretty(B);
            fmpz_mod_mat_print_pretty(C);
            fmpz_mod_mat_print_pretty(D);
            abort();
        }

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
        fmpz_mod_mat_clear(D);

        fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
#else
    FLINT_TEST_CLEANUP(state);

    flint_printf("SKIPPED\n");
#endif
    return 0;
}

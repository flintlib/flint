/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_mul_classical_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    slong i, max_threads = 5;
    slong tmul = 1000;
#ifdef _WIN32
    tmul = 50;
#endif

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, C, D;
        fmpz_mod_ctx_t ctx;

        slong m, k, n;

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        fmpz_mod_ctx_init_rand_bits(ctx, state, 100);

        fmpz_mod_mat_init(A, m, n, ctx);
        fmpz_mod_mat_init(B, n, k, ctx);
        fmpz_mod_mat_init(C, m, k, ctx);
        fmpz_mod_mat_init(D, m, k, ctx);

        fmpz_mod_mat_randtest(A, state, ctx);
        fmpz_mod_mat_randtest(B, state, ctx);
        fmpz_mod_mat_randtest(C, state, ctx);  /* make sure noise in the output is ok */

        fmpz_mod_mat_mul_classical_threaded(C, A, B, ctx);
        fmpz_mod_mat_mul(D, A, B, ctx);

        if (!fmpz_mod_mat_equal(C, D, ctx))
        {
            flint_printf("FAIL: results not equal\n");
            fmpz_mod_mat_print_pretty(A, ctx);
            fmpz_mod_mat_print_pretty(B, ctx);
            fmpz_mod_mat_print_pretty(C, ctx);
            fmpz_mod_mat_print_pretty(D, ctx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A, ctx);
        fmpz_mod_mat_clear(B, ctx);
        fmpz_mod_mat_clear(C, ctx);
        fmpz_mod_mat_clear(D, ctx);

        fmpz_mod_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}

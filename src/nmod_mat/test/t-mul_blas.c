/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

/* generate a worst case matrix for blas */
void nmod_mat_randfull_half(nmod_mat_t mat, flint_rand_t state)
{
    slong i, j;
    slong r = mat->r;
    slong c = mat->c;

    for (i = 0; i < r; i++)
    for (j = 0; j < c; j++)
    {
        mat->rows[i][j] = mat->mod.n/2;
        if (mat->mod.n > 2 && (mat->mod.n % 2))
            mat->rows[i][j] += n_randint(state, 2);
    }
}

TEST_FUNCTION_START(nmod_mat_mul_blas, state)
{
    slong i, max_threads = 5;

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t modulus;
        slong m, k, n;

        m = n_randint(state, 150) + 2;
        k = n_randint(state, 150) + 2;
        n = n_randint(state, 150) + 2;

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(state, 3))
        {
            case 0:
                modulus = n_randtest_not_zero(state);
                break;
            case 1:
                modulus = UWORD_MAX/2 + 1 - n_randbits(state, 4);
                break;
            default:
                modulus = UWORD_MAX - n_randbits(state, 4);
                break;
        }

        nmod_mat_init(A, m, n, modulus);
        nmod_mat_init(B, n, k, modulus);
        nmod_mat_init(C, m, k, modulus);
        nmod_mat_init(D, m, k, modulus);

        if (n_randint(state, 2))
            nmod_mat_randfull_half(A, state);
        else
            nmod_mat_randfull(A, state);

        if (n_randint(state, 2))
            nmod_mat_randfull_half(B, state);
        else
            nmod_mat_randfull(B, state);

        nmod_mat_randtest(C, state);

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        if (nmod_mat_mul_blas(C, A, B))
        {
            nmod_mat_mul_classical(D, A, B);

            if (!nmod_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal\n");
                flint_printf("m: %wd, k: %wd, n: %wd, mod: %wu\n", m, k, n, modulus);
                fflush(stdout);
                flint_abort();
            }
        }
#if FLINT_USES_BLAS && FLINT_BITS == 64
        else
        {
            flint_printf("FAIL: blas should have worked\n");
            fflush(stdout);
            flint_abort();
        }
#endif

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}

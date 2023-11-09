/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_mat.h"

/* Defined in t-mul.c and t-mul_classical_threaded.c */
#ifndef nmod_mat_mul_check
#define nmod_mat_mul_check nmod_mat_mul_check
void
nmod_mat_mul_check(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, k;

    mp_limb_t s0, s1, s2;
    mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s0 = s1 = s2 = UWORD(0);

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, A->rows[i][k], B->rows[k][j]);
                add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
            }

            NMOD_RED(s2, s2, C->mod);
            NMOD_RED3(s0, s2, s1, s0, C->mod);
            C->rows[i][j] = s0;
        }
    }
}
#endif

TEST_FUNCTION_START(nmod_mat_mul, state)
{
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t mod;

        slong m, k, n;

        m = n_randint(state, 75);
        k = n_randint(state, 75);
        n = n_randint(state, 75);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(state, 3))
        {
            case 0:
                mod = n_randtest_not_zero(state);
                break;
            case 1:
                mod = UWORD_MAX/2 + 1 - n_randbits(state, 4);
                break;
            case 2:
            default:
                mod = UWORD_MAX - n_randbits(state, 4);
                break;
        }

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);

        if (n_randint(state, 2))
            nmod_mat_randtest(A, state);
        else
            nmod_mat_randfull(A, state);

        if (n_randint(state, 2))
            nmod_mat_randtest(B, state);
        else
            nmod_mat_randfull(B, state);

        nmod_mat_randtest(C, state);  /* make sure noise in the output is ok */

        nmod_mat_mul(C, A, B);
        nmod_mat_mul_check(D, A, B);

        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            fflush(stdout);
            flint_abort();
        }

        if (n == k)
        {
            nmod_mat_mul(A, A, B);

            if (!nmod_mat_equal(A, C))
            {
                flint_printf("FAIL: aliasing failed\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    /* Test aliasing with windows */
    {
        nmod_mat_t A, B, A_window;

        nmod_mat_init(A, 2, 2, 3);
        nmod_mat_init(B, 2, 2, 3);

        nmod_mat_window_init(A_window, A, 0, 0, 2, 2);

        nmod_mat_one(A);
        nmod_mat_one(B);
        nmod_mat_entry(B, 0, 1) = 1;
        nmod_mat_entry(B, 1, 0) = 1;

        nmod_mat_mul(A_window, B, A_window);

        if (!nmod_mat_equal(A, B))
        {
            flint_printf("FAIL: window aliasing failed\n");
            nmod_mat_print_pretty(A); flint_printf("\n\n");
            nmod_mat_print_pretty(B); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_window_clear(A_window);
        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}

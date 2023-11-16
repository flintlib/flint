/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "bool_mat.h"
#include "perm.h"

/* Defined in t-all_pairs_longest_walk.c and t-transitive_closure.c */
#ifndef _bool_mat_permute
#define _bool_mat_permute _bool_mat_permute
/* permute rows and columns of a square matrix */
void
_bool_mat_permute(bool_mat_t B, const bool_mat_t A, const slong *perm)
{
    slong n, i, j;
    if ((A == B) || !bool_mat_is_square(A)) flint_abort(); /* assert */
    n = bool_mat_nrows(A);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            bool_mat_set_entry(
                    B, perm[i], perm[j], bool_mat_get_entry(A, i, j));
        }
    }
}
#endif

void
_bool_mat_transitive_closure_powering(bool_mat_t B, const bool_mat_t A)
{
    slong n, k;
    bool_mat_t T;
    if ((A == B) || !bool_mat_is_square(A)) flint_abort(); /* assert */
    n = bool_mat_nrows(A);
    bool_mat_init(T, n, n);
    bool_mat_one(T);
    bool_mat_zero(B);
    for (k = 0; k < n; k++)
    {
        bool_mat_mul(T, T, A);
        bool_mat_add(B, B, T);
    }
    bool_mat_clear(T);
}

TEST_FUNCTION_START(bool_mat_transitive_closure, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m;
        bool_mat_t A, B, C, D;

        m = n_randint(state, 50);

        bool_mat_init(A, m, m);
        bool_mat_init(B, m, m);
        bool_mat_init(C, m, m);
        bool_mat_init(D, m, m);

        bool_mat_randtest(A, state);
        bool_mat_randtest(B, state);

        bool_mat_transitive_closure(B, A);

        /* transitivity */
        if (!bool_mat_is_transitive(B))
        {
            flint_printf("FAIL (transitivity)\n");
            bool_mat_print(A); flint_printf("\n");
            bool_mat_print(B); flint_printf("\n");
            flint_abort();
        }

        /* monotonicity */
        {
            bool_mat_complement(C, B);
            bool_mat_mul_entrywise(D, A, C);
            if (bool_mat_any(D))
            {
                flint_printf("FAIL (monotonicity)\n");
                bool_mat_print(A); flint_printf("\n");
                bool_mat_print(B); flint_printf("\n");
                flint_abort();
            }
        }

        /* aliasing */
        {
            bool_mat_set(C, A);
            bool_mat_transitive_closure(C, C);
            if (!bool_mat_equal(B, C))
            {
                flint_printf("FAIL (aliasing)\n");
                bool_mat_print(A); flint_printf("\n");
                bool_mat_print(B); flint_printf("\n");
                bool_mat_print(C); flint_printf("\n");
                flint_abort();
            }
        }

        /* test commutativity of permutation with transitive closure */
        {
            slong *perm;
            perm = flint_malloc(m * sizeof(slong));
            _perm_randtest(perm, m, state);

            /* C is the transitive closure of the permutation of A */
            bool_mat_randtest(C, state);
            _bool_mat_permute(C, A, perm);
            bool_mat_transitive_closure(C, C);

            /* D is the permutation of the transitive closure of A */
            bool_mat_randtest(D, state);
            _bool_mat_permute(D, B, perm);

            if (!bool_mat_equal(C, D))
            {
                flint_printf("FAIL (commutativity with permutation)\n");
                bool_mat_print(A); flint_printf("\n");
                bool_mat_print(B); flint_printf("\n");
                bool_mat_print(C); flint_printf("\n");
                bool_mat_print(D); flint_printf("\n");
                flint_abort();
            }
            flint_free(perm);
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
        bool_mat_clear(C);
        bool_mat_clear(D);
    }

    /* check transitive closure using brute force with smallish matrices */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m;
        bool_mat_t A, B, C;

        m = n_randint(state, 10);

        bool_mat_init(A, m, m);
        bool_mat_init(B, m, m);
        bool_mat_init(C, m, m);

        bool_mat_randtest(A, state);
        bool_mat_randtest(B, state);
        bool_mat_randtest(C, state);

        bool_mat_transitive_closure(B, A);
        _bool_mat_transitive_closure_powering(C, A);

        if (!bool_mat_equal(B, C))
        {
            flint_printf("FAIL (brute force)\n");
            bool_mat_print(A); flint_printf("\n");
            bool_mat_print(B); flint_printf("\n");
            bool_mat_print(C); flint_printf("\n");
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
        bool_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}

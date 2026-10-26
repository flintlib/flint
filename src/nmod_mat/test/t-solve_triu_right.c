/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

/* Defined in t-solve_tril_right.c and t-solve_triu_right.c */
#ifndef _solve_tri_right_check
#define _solve_tri_right_check _solve_tri_right_check
/*
    Checks the solution of X T = B for a random r x r triangular T (lower
    if lower is nonzero) and k x r matrix X, with the algorithm algo
    (0: default, 1: classical, 2: recursive), in windows, with aliasing,
    and with garbage in the entries of T that must not be read (the other
    triangular part, and the diagonal if unit). Returns 0 if correct.
*/
static int
_solve_tri_right_check(flint_rand_t state, int lower, int algo)
{
    nmod_mat_t T, X, B, Y, BB, YY, W, WB;
    slong r, k, i, j;
    ulong p;
    int unit, res = 0;

    p = n_randtest_prime(state, 0);
    r = n_randint(state, 100);
    k = n_randint(state, 100);
    if (n_randint(state, 4) == 0)
    {
        r = n_randint(state, 200);
        k = n_randint(state, 200);
    }
    unit = n_randint(state, 2);

    nmod_mat_init(T, r, r, p);
    nmod_mat_init(X, k, r, p);
    nmod_mat_init(B, k, r, p);
    nmod_mat_init(Y, k, r, p);

    if (lower)
        nmod_mat_randtril(T, state, unit);
    else
        nmod_mat_randtriu(T, state, unit);
    nmod_mat_randtest(X, state);
    nmod_mat_mul(B, X, T);

    /* entries that must not be read */
    for (i = 0; i < r; i++)
        for (j = 0; j < r; j++)
            if ((lower ? j > i : j < i) || (unit && i == j))
                nmod_mat_entry(T, i, j) = n_randint(state, p);

#define SOLVE(X, T, B) \
    do { \
        if (lower) \
        { \
            if (algo == 0) nmod_mat_solve_tril_right(X, T, B, unit); \
            else if (algo == 1) nmod_mat_solve_tril_right_classical(X, T, B, unit); \
            else nmod_mat_solve_tril_right_recursive(X, T, B, unit); \
        } \
        else \
        { \
            if (algo == 0) nmod_mat_solve_triu_right(X, T, B, unit); \
            else if (algo == 1) nmod_mat_solve_triu_right_classical(X, T, B, unit); \
            else nmod_mat_solve_triu_right_recursive(X, T, B, unit); \
        } \
    } while (0)

    SOLVE(Y, T, B);
    if (!nmod_mat_equal(Y, X))
        res = 1;

    /* aliasing */
    nmod_mat_set(Y, B);
    SOLVE(Y, T, Y);
    if (!nmod_mat_equal(Y, X))
        res = 2;

    /* windows of larger matrices */
    nmod_mat_init(BB, k + 2, r + 3, p);
    nmod_mat_init(YY, k + 1, r + 2, p);
    nmod_mat_randtest(BB, state);
    nmod_mat_randtest(YY, state);
    nmod_mat_window_init(WB, BB, 1, 2, k + 1, r + 2);
    nmod_mat_window_init(W, YY, 0, 1, k, r + 1);
    nmod_mat_set(WB, B);
    SOLVE(W, T, WB);
    if (!nmod_mat_equal(W, X))
        res = 3;
    nmod_mat_window_clear(W);
    nmod_mat_window_clear(WB);

#undef SOLVE

    nmod_mat_clear(T);
    nmod_mat_clear(X);
    nmod_mat_clear(B);
    nmod_mat_clear(Y);
    nmod_mat_clear(BB);
    nmod_mat_clear(YY);

    return res;
}
#endif

TEST_FUNCTION_START(nmod_mat_solve_triu_right, state)
{
    slong iter;
    int algo, res;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        for (algo = 0; algo < 3; algo++)
        {
            res = _solve_tri_right_check(state, 0, algo);
            if (res != 0)
                TEST_FUNCTION_FAIL("algorithm %d, error %d\n", algo, res);
        }
    }

    TEST_FUNCTION_END(state);
}

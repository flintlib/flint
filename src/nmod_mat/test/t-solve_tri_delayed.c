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
#include "nmod_mat/impl.h"

/* random modulus, with moduli close to 2^(FLINT_BITS / 2) (one-limb
   accumulators reduced along the way), 2^(FLINT_BITS - 1) and
   2^FLINT_BITS (two and three limbs) */
static ulong
_tri_modulus(flint_rand_t state)
{
    switch (n_randint(state, 6))
    {
        case 0:
            return n_randtest_not_zero(state);
        case 1:
            return (UWORD(1) << (FLINT_BITS / 2)) - n_randint(state, 1000);
        case 2:
            return (UWORD(1) << (FLINT_BITS / 2 - 2)) + n_randint(state, 1000);
        case 3:
            return (UWORD(1) << (FLINT_BITS - 1)) + n_randint(state, 1000);
        case 4:
            return UWORD_MAX - n_randint(state, 1000);
        default:
            return n_randtest_prime(state, 0);
    }
}

/*
    X T = B (right) or T X = B with the algorithm algo (0: rows, 1: axpy,
    right only, 2: classical), for a random r x r triangular T with
    invertible diagonal, garbage in the entries that must not be read;
    B = X0 T (or T X0), or arbitrary B and the solution checked by
    multiplication, in windows, with aliasing. Returns 0 if correct.
*/
static int
_check(flint_rand_t state, int algo)
{
    nmod_mat_t T, Tc, X0, B, Y, Z, BB, YY, WB, WY;
    slong r, k, i, j, br, bc;
    ulong n;
    int unit, upper, right, byproduct, res = 0;

    n = _tri_modulus(state);
    r = n_randint(state, 40);
    k = n_randint(state, 40);
    if (n_randint(state, 8) == 0)
        r = n_randint(state, 150);
    if (n_randint(state, 8) == 0)
        k = n_randint(state, 600);
    unit = n_randint(state, 2);
    upper = n_randint(state, 2);
    right = (algo == 1) ? 1 : n_randint(state, 2);
    byproduct = n_randint(state, 2);

    br = right ? k : r;
    bc = right ? r : k;

    /* Tc: the triangular matrix, T: with garbage */
    nmod_mat_init(Tc, r, r, n);
    nmod_mat_init(T, r, r, n);
    nmod_mat_init(X0, br, bc, n);
    nmod_mat_init(B, br, bc, n);
    nmod_mat_init(Y, br, bc, n);
    nmod_mat_init(Z, br, bc, n);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            ulong e, g;

            if (i == j)
            {
                do
                    e = n_randint(state, n);
                while (n_gcd(e, n) != 1);
                nmod_mat_entry(Tc, i, i) = unit ? 1 : e;
                nmod_mat_entry(T, i, i) = unit ? n_randint(state, n) : e;
            }
            else if (upper ? (j > i) : (j < i))
            {
                e = (n_randint(state, 4) == 0) ? n - 1 : n_randint(state, n);
                nmod_mat_entry(Tc, i, j) = e;
                nmod_mat_entry(T, i, j) = e;
            }
            else
            {
                g = n_randint(state, n);
                nmod_mat_entry(T, i, j) = g;
            }
        }
    }

    if (byproduct)
    {
        nmod_mat_randtest(X0, state);
        if (right)
            nmod_mat_mul(B, X0, Tc);
        else
            nmod_mat_mul(B, Tc, X0);
    }
    else
    {
        if (n_randint(state, 2))
            nmod_mat_randtest(B, state);
        else
            for (i = 0; i < br; i++)
                for (j = 0; j < bc; j++)
                    nmod_mat_entry(B, i, j) = n - 1;
    }

#define SOLVE(Y, B) \
    do { \
        if (algo == 0) \
            _nmod_mat_solve_tri_rows_delayed(Y, T, B, unit, upper, right); \
        else if (algo == 1) \
            _nmod_mat_solve_tri_right_axpy_delayed(Y, T, B, unit, upper); \
        else if (right && upper) \
            nmod_mat_solve_triu_right_classical(Y, T, B, unit); \
        else if (right) \
            nmod_mat_solve_tril_right_classical(Y, T, B, unit); \
        else if (upper) \
            nmod_mat_solve_triu_classical(Y, T, B, unit); \
        else \
            nmod_mat_solve_tril_classical(Y, T, B, unit); \
    } while (0)

#define CHECK(Y, code) \
    do { \
        if (byproduct) \
        { \
            if (!nmod_mat_equal(Y, X0)) \
                res = code; \
        } \
        else \
        { \
            if (right) \
                nmod_mat_mul(Z, Y, Tc); \
            else \
                nmod_mat_mul(Z, Tc, Y); \
            if (!nmod_mat_equal(Z, B)) \
                res = code; \
        } \
    } while (0)

    SOLVE(Y, B);
    CHECK(Y, 1);

    /* aliasing */
    nmod_mat_set(Y, B);
    SOLVE(Y, Y);
    CHECK(Y, 2);

    /* windows */
    nmod_mat_init(BB, br + 2, bc + 3, n);
    nmod_mat_init(YY, br + 3, bc + 1, n);
    nmod_mat_randtest(BB, state);
    nmod_mat_randtest(YY, state);
    nmod_mat_window_init(WB, BB, 1, 2, br + 1, bc + 2);
    nmod_mat_window_init(WY, YY, 2, 1, br + 2, bc + 1);
    nmod_mat_set(WB, B);
    SOLVE(WY, WB);
    nmod_mat_set(Y, WY);
    CHECK(Y, 3);
    nmod_mat_window_clear(WB);
    nmod_mat_window_clear(WY);

#undef SOLVE
#undef CHECK

    nmod_mat_clear(T);
    nmod_mat_clear(Tc);
    nmod_mat_clear(X0);
    nmod_mat_clear(B);
    nmod_mat_clear(Y);
    nmod_mat_clear(Z);
    nmod_mat_clear(BB);
    nmod_mat_clear(YY);

    return res;
}

TEST_FUNCTION_START(nmod_mat_solve_tri_delayed, state)
{
    slong iter;
    int algo, res;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        for (algo = 0; algo < 3; algo++)
        {
            res = _check(state, algo);
            if (res != 0)
                TEST_FUNCTION_FAIL("algorithm %d, error %d\n", algo, res);
        }
    }

    TEST_FUNCTION_END(state);
}

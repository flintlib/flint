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
#include "perm.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"
#include "rank_profile_helpers.h"

/* cutoff of the recursion in the tests, chosen at random */
static slong _ple_test_cutoff = 1;

static slong
_ple_recursive_test(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    return _nmod_mat_ple_recursive_cutoff(P, Q, A, max_nonpivots, _ple_test_cutoff);
}

/* the algorithms, called through their public interface */
#define PLE_NUM_ALGOS 6

static slong
_ple_algo(int algo, slong * P, slong * Q, nmod_mat_t A, int flags)
{
    switch (algo)
    {
        case 0:
            return nmod_mat_ple_classical(P, Q, A, flags);
        case 1:
            return nmod_mat_ple_classical_delayed(P, Q, A, flags);
        case 2:
            return nmod_mat_ple_crout(P, Q, A, flags);
        case 3:
            return nmod_mat_ple_recursive(P, Q, A, flags);
        case 4:
            return _nmod_mat_ple_wrap(P, Q, A, flags, _ple_recursive_test);
        default:
            return nmod_mat_ple(P, Q, A, flags);
    }
}

static const char * _ple_algo_name[PLE_NUM_ALGOS] =
    { "classical", "classical_delayed", "crout", "recursive",
      "recursive (small cutoff)", "default" };

/* the internal functions, with a bound max_nonpivots */
#define PLE_NUM_INTERNAL 5

static slong
_ple_internal(int algo, slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots)
{
    switch (algo)
    {
        case 0:
            return _nmod_mat_ple_classical(P, Q, A, max_nonpivots);
        case 1:
            return _nmod_mat_ple_classical_delayed(P, Q, A, max_nonpivots);
        case 2:
            return _nmod_mat_ple_crout(P, Q, A, max_nonpivots);
        case 3:
            return _nmod_mat_ple_recursive(P, Q, A, max_nonpivots);
        default:
            return _ple_recursive_test(P, Q, A, max_nonpivots);
    }
}

/* the number of columns without pivot met while some rows are not pivot
   rows yet, from the canonical output */
static slong
_ple_nonpivots_met(const slong * Q, slong m, slong n, slong rank)
{
    if (rank < m)
        return n - rank;
    if (m == 0)
        return 0;
    return Q[m - 1] - (m - 1);
}

TEST_FUNCTION_START(nmod_mat_ple, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_mat_t A, LU, LU0;
        slong m, n, rank, exp_rank, r;
        slong * P, * Q, * P0, * Q0, * mrp_r, * mrp_c;
        int kind, algo, res, known, check_rpm;
        ulong p;

        m = n_randint(state, 13);
        n = n_randint(state, 13);
        if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 60);
            n = n_randint(state, 60);
        }
        if (n_randint(state, 100) == 0)
        {
            m = n_randint(state, 200);
            n = n_randint(state, 200);
        }
        p = _rp_rand_prime(state);
        check_rpm = (m <= 12 && n <= 12);
        _ple_test_cutoff = 1 + n_randint(state, 10);

        nmod_mat_init(A, m, n, p);
        mrp_r = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        mrp_c = flint_malloc((FLINT_MIN(m, n) + 1) * sizeof(slong));
        exp_rank = -1;
        known = 0;

        kind = n_randint(state, 5);
        if (kind == 0 || kind == 1)
        {
            /* prescribed rank profile matrix; kind 1: low rank */
            exp_rank = n_randint(state, FLINT_MIN(m, n) + 1);
            if (kind == 1)
                exp_rank = n_randint(state, FLINT_MIN(exp_rank, 3) + 1);
            _rp_rand_positions(mrp_r, mrp_c, m, n, exp_rank, state);
            nmod_mat_randrank_profile(A, state, mrp_r, mrp_c, exp_rank);
            known = 1;
        }
        else if (kind == 2)
        {
            nmod_mat_randtest(A, state);
        }
        else if (kind == 3)
        {
            /* largest growth of the unreduced entries: A = L U with L unit
               lower triangular with ones below the diagonal and U upper
               triangular with all entries p - 1 */
            slong i, j;

            for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                    nmod_mat_entry(A, i, j) = nmod_mul(
                            nmod_set_ui(FLINT_MIN(i, j) + 1, A->mod), p - 1, A->mod);
            exp_rank = FLINT_MIN(m, n);
            for (i = 0; i < exp_rank; i++)
                mrp_r[i] = mrp_c[i] = i;
            known = 1;
        }
        else
        {
            r = n_randint(state, FLINT_MIN(m, n) + 1);
            nmod_mat_randrank(A, state, r);
            if (n_randint(state, 2))
                nmod_mat_randops(A, state, n_randint(state, 2 * m * n + 1));
            exp_rank = r;
        }

        P = flint_malloc((m + 1) * sizeof(slong));
        P0 = flint_malloc((m + 1) * sizeof(slong));
        Q = flint_malloc((n + 1) * sizeof(slong));
        Q0 = flint_malloc((n + 1) * sizeof(slong));
        nmod_mat_init(LU, m, n, p);
        nmod_mat_init(LU0, m, n, p);

        /* each algorithm: canonical output, identical for all */
        for (algo = 0; algo < PLE_NUM_ALGOS; algo++)
        {
            nmod_mat_set(LU, A);
            rank = _ple_algo(algo, P, Q, LU, 0);

            res = _rp_check_ple(P, Q, LU, rank, A, known ? mrp_r : NULL,
                                mrp_c, exp_rank, check_rpm);
            if (res == 0 && exp_rank >= 0 && rank != exp_rank)
                res = 100;

            if (res != 0)
                TEST_FUNCTION_FAIL(
                        "%s: error %d\n"
                        "m = %wd, n = %wd, p = %wu, rank = %wd, expected %wd\n"
                        "A = %{nmod_mat}\nLU = %{nmod_mat}\n"
                        "P = %{slong*}\nQ = %{slong*}\n",
                        _ple_algo_name[algo], res, m, n, p, rank, exp_rank,
                        A, LU, P, m, Q, n);

            if (algo == 0)
            {
                nmod_mat_set(LU0, LU);
                _perm_set(P0, P, m);
                _perm_set(Q0, Q, n);
                exp_rank = rank;
            }
            else if (!nmod_mat_equal(LU, LU0) || !_perm_equal(P, P0, m)
                        || !_perm_equal(Q, Q0, n))
                TEST_FUNCTION_FAIL(
                        "%s: output differs from classical\n"
                        "A = %{nmod_mat}\nLU = %{nmod_mat}\nLU0 = %{nmod_mat}\n",
                        _ple_algo_name[algo], A, LU, LU0);

            /* rank check: 0 if the rank is not min(m, n), else the same */
            nmod_mat_set(LU, A);
            r = _ple_algo(algo, P, Q, LU, NMOD_MAT_PLE_RANK_CHECK);
            if (rank < FLINT_MIN(m, n) ? r != 0
                    : (r != rank || !nmod_mat_equal(LU, LU0)
                        || !_perm_equal(P, P0, m) || !_perm_equal(Q, Q0, n)))
                TEST_FUNCTION_FAIL(
                        "%s: rank check, returned %wd, rank %wd\n"
                        "A = %{nmod_mat}\n",
                        _ple_algo_name[algo], r, rank, A);

            /* pivots_nonpivots NULL */
            nmod_mat_set(LU, A);
            r = _ple_algo(algo, P, NULL, LU, 0);
            if (r != rank || !nmod_mat_equal(LU, LU0) || !_perm_equal(P, P0, m))
                TEST_FUNCTION_FAIL("%s: NULL pivots_nonpivots\n",
                                   _ple_algo_name[algo]);

            /* in a window: the rest of the matrix is not modified */
            {
                nmod_mat_t B, B0, W;
                slong i, j;

                nmod_mat_init(B, m + 3, n + 2, p);
                nmod_mat_randtest(B, state);
                nmod_mat_window_init(W, B, 2, 1, m + 2, n + 1);
                nmod_mat_set(W, A);
                nmod_mat_init_set(B0, B);
                r = _ple_algo(algo, P, Q, W, 0);
                if (r != rank || !nmod_mat_equal(W, LU0)
                        || !_perm_equal(P, P0, m) || !_perm_equal(Q, Q0, n))
                    TEST_FUNCTION_FAIL("%s: window\n", _ple_algo_name[algo]);
                for (i = 0; i < m + 3; i++)
                    for (j = 0; j < n + 2; j++)
                        if ((i < 2 || i >= m + 2 || j < 1 || j >= n + 1)
                                && nmod_mat_entry(B, i, j) != nmod_mat_entry(B0, i, j))
                            TEST_FUNCTION_FAIL("%s: window, entry (%wd, %wd) modified\n",
                                               _ple_algo_name[algo], i, j);
                nmod_mat_window_clear(W);
                nmod_mat_clear(B);
                nmod_mat_clear(B0);
            }
        }

        /* internal functions: -1 exactly when more than max_nonpivots
           columns without pivot are met, otherwise the same output */
        for (algo = 0; algo < PLE_NUM_INTERNAL; algo++)
        {
            slong met = _ple_nonpivots_met(Q0, m, n, exp_rank);
            slong bound = (slong) n_randint(state, met + 3) - 1;

            nmod_mat_set(LU, A);
            r = _ple_internal(algo, P, Q, LU, bound);
            if ((bound >= 0 && met > bound) ? r != -1
                    : (r != exp_rank || !nmod_mat_equal(LU, LU0)
                        || !_perm_equal(P, P0, m) || !_perm_equal(Q, Q0, n)))
                TEST_FUNCTION_FAIL(
                        "internal %d: bound %wd, %wd met, returned %wd, rank %wd\n"
                        "A = %{nmod_mat}\n",
                        algo, bound, met, r, exp_rank, A);
        }

        /* same rank and pivot columns as nmod_mat_lu */
        {
            slong * Plu = flint_malloc((m + 1) * sizeof(slong));
            slong * Qlu = flint_malloc((n + 1) * sizeof(slong));
            nmod_mat_t LU2;

            nmod_mat_init_set(LU2, A);
            r = nmod_mat_lu_with_pivots(Plu, Qlu, LU2);
            if (r != exp_rank || !_perm_equal(Qlu, Q0, n))
                TEST_FUNCTION_FAIL("rank or pivots differ from nmod_mat_lu\n");
            nmod_mat_clear(LU2);
            flint_free(Plu);
            flint_free(Qlu);
        }

        nmod_mat_clear(A);
        nmod_mat_clear(LU);
        nmod_mat_clear(LU0);
        flint_free(P);
        flint_free(Q);
        flint_free(P0);
        flint_free(Q0);
        flint_free(mrp_r);
        flint_free(mrp_c);
    }

    TEST_FUNCTION_END(state);
}

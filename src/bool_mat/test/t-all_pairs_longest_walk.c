/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "bool_mat.h"

int
_is_superficially_ok_entrywise(const fmpz_mat_t A)
{
    slong n, i, j, d;

    if (!fmpz_mat_is_square(A))
        return 0;

    n = fmpz_mat_nrows(A);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            d = fmpz_get_si(fmpz_mat_entry(A, i, j));

            if (d < -2 || d > n-1)
                return 0;

            if (i == j && d == -1)
                return 0;

            if (i != j && d == 0)
                return 0;
        }
    }
    return 1;
}

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
_fmpz_mat_permute(fmpz_mat_t B, const fmpz_mat_t A, const slong *perm)
{
    slong n, i, j;
    if (!fmpz_mat_is_square(A)) flint_abort(); /* assert */
    if (A == B) flint_abort(); /* assert */
    n = fmpz_mat_nrows(A);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_set(
                    fmpz_mat_entry(B, perm[i], perm[j]),
                    fmpz_mat_entry(A, i, j));
        }
    }
}

void
_brute_force_all_pairs_longest_walk(fmpz_mat_t B, const bool_mat_t A)
{
    slong i, j, n;

    n = bool_mat_nrows(A);

    /* set entries of B according to the longest observed walk */
    {
        slong k;
        bool_mat_t T;
        bool_mat_init(T, n, n);
        bool_mat_one(T);
        fmpz_mat_zero(B);
        for (k = 0; k < 2*n+1; k++)
        {
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (bool_mat_get_entry(T, i, j))
                    {
                        fmpz_set_si(fmpz_mat_entry(B, i, j), k);
                    }
                }
            }
            bool_mat_mul(T, T, A);
        }
        bool_mat_clear(T);
    }

    /* set special values 0, -1, -2 */
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                slong x;
                fmpz *p;
                p = fmpz_mat_entry(B, i, j);
                x = fmpz_get_si(p);
                if (x < 1)
                {
                    x = (i == j) ? 0 : -1;
                }
                else if (x > n-1)
                {
                    x = -2;
                }
                fmpz_set_si(p, x);
            }
        }
    }
}

TEST_FUNCTION_START(bool_mat_all_pairs_longest_walk, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, degree;
        bool_mat_t A;
        fmpz_mat_t B, C, D;

        m = n_randint(state, 50);

        bool_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, m, m);

        bool_mat_randtest(A, state);

        degree = bool_mat_all_pairs_longest_walk(B, A);

        /* entrywise reasonability */
        if (!_is_superficially_ok_entrywise(B))
        {
            flint_printf("FAIL (entrywise)\n");
            bool_mat_print(A); flint_printf("\n");
            fmpz_mat_print_pretty(B); flint_printf("\n");
            flint_abort();
        }

        /* nilpotency degree */
        {
            slong nildegree = bool_mat_nilpotency_degree(A);
            if (nildegree != degree + 1)
            {
                flint_printf("FAIL (nilpotency degree)\n");
                bool_mat_print(A); flint_printf("\n");
                fmpz_mat_print_pretty(B); flint_printf("\n");
                flint_printf("nildegree=%wd degree=%wd\n", nildegree, degree);
                flint_abort();
            }
        }

        /* transitive closure */
        {
            slong i, j;
            bool_mat_t U, V;

            bool_mat_init(U, m, m);
            bool_mat_transitive_closure(U, A);

            bool_mat_init(V, m, m);
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    slong x = fmpz_get_si(fmpz_mat_entry(B, i, j));
                    bool_mat_set_entry(V, i, j, (x != -1 && x != 0));
                }
            }

            if (!bool_mat_equal(U, V))
            {
                flint_printf("FAIL (transitive closure)\n");
                bool_mat_print(A); flint_printf("\n");
                fmpz_mat_print_pretty(B); flint_printf("\n");
                bool_mat_print(U); flint_printf("\n");
                bool_mat_print(V); flint_printf("\n");
                flint_abort();
            }

            bool_mat_clear(U);
            bool_mat_clear(V);
        }

        /* test commutativity of all-pairs-longest-walk with permutation */
        {
            bool_mat_t Ap;
            slong *perm;

            bool_mat_init(Ap, m, m);
            perm = flint_malloc(m * sizeof(slong));
            _perm_randtest(perm, m, state);

            /* C is the all-pairs-longest-walk of the permutation of A */
            _bool_mat_permute(Ap, A, perm);
            bool_mat_all_pairs_longest_walk(C, Ap);

            /* D is the permutation of the all-pairs-longest-walk of A */
            _fmpz_mat_permute(D, B, perm);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL (permutation)\n");
                bool_mat_print(A); flint_printf("\n");
                fmpz_mat_print_pretty(B); flint_printf("\n");
                fmpz_mat_print_pretty(C); flint_printf("\n");
                fmpz_mat_print_pretty(D); flint_printf("\n");
                flint_abort();
            }

            flint_free(perm);
            bool_mat_clear(Ap);
        }

        bool_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    /* use powering to check small random examples */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m;
        bool_mat_t A;
        fmpz_mat_t B, C;

        m = n_randint(state, 10);

        bool_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);

        bool_mat_randtest(A, state);
        bool_mat_all_pairs_longest_walk(B, A);

        _brute_force_all_pairs_longest_walk(C, A);

        if (!fmpz_mat_equal(B, C))
        {
            flint_printf("FAIL (powering)\n");
            bool_mat_print(A); flint_printf("\n");
            fmpz_mat_print_pretty(B); flint_printf("\n");
            fmpz_mat_print_pretty(C); flint_printf("\n");
            flint_abort();
        }

        bool_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    /* special matrices */
    {
        slong m;
        for (m = 1; m < 100; m++)
        {
            slong i, j, length;
            bool_mat_t A;
            fmpz_mat_t B;

            bool_mat_init(A, m, m);
            fmpz_mat_init(B, m, m);

            /* directed path */
            {
                bool_mat_directed_path(A);

                length = bool_mat_all_pairs_longest_walk(B, A);
                if (length != m-1)
                {
                    flint_printf("FAIL (directed path)\n");
                    bool_mat_print(A); flint_printf("\n");
                    fmpz_mat_print_pretty(B); flint_printf("\n");
                    flint_printf("m=%wd length=%wd\n", m, length);
                    flint_abort();
                }

                for (i = 0; i < m; i++)
                {
                    for (j = 0; j < m; j++)
                    {
                        if (fmpz_get_si(fmpz_mat_entry(B, i, j)) !=
                            FLINT_MAX(-1, j - i))
                        {
                            flint_printf("FAIL (directed path)\n");
                            bool_mat_print(A); flint_printf("\n");
                            fmpz_mat_print_pretty(B); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }

            /* directed cycle */
            {
                bool_mat_directed_cycle(A);

                length = bool_mat_all_pairs_longest_walk(B, A);
                if (length != -2)
                {
                    flint_printf("FAIL (directed cycle)\n");
                    bool_mat_print(A); flint_printf("\n");
                    fmpz_mat_print_pretty(B); flint_printf("\n");
                    flint_printf("m=%wd length=%wd\n", m, length);
                    flint_abort();
                }

                for (i = 0; i < m; i++)
                {
                    for (j = 0; j < m; j++)
                    {
                        if (fmpz_get_si(fmpz_mat_entry(B, i, j)) != -2)
                        {
                            flint_printf("FAIL (directed cycle)\n");
                            bool_mat_print(A); flint_printf("\n");
                            fmpz_mat_print_pretty(B); flint_printf("\n");
                            flint_abort();
                        }
                    }
                }
            }

            bool_mat_clear(A);
            fmpz_mat_clear(B);
        }
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

/*
    Helper: build m x n matrix from slong array, compute snf_transform,
    check S matches expected diagonal (verified with Magma), and verify
    U*A*V == S with |det(U)| == |det(V)| == 1.

    Expected SNF entries are passed as decimal strings so test cases
    with entries exceeding slong range are portable across 32/64-bit.
*/
static void
_test_transform_case(slong m, slong n, const slong * data,
    const char * const * expected_snf)
{
    fmpz_mat_t A, S, U, V, T1, T2;
    fmpz_t det, expected;
    slong i, j, d;

    d = FLINT_MIN(m, n);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(S, m, n);
    fmpz_mat_init(U, m, m);
    fmpz_mat_init(V, n, n);
    fmpz_mat_init(T1, m, n);
    fmpz_mat_init(T2, m, n);
    fmpz_init(det);
    fmpz_init(expected);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), data[i * n + j]);

    fmpz_mat_snf_transform(S, U, V, A);

    /* Check diagonal matches expected */
    for (i = 0; i < d; i++)
    {
        fmpz_set_str(expected, expected_snf[i], 10);
        if (!fmpz_equal(fmpz_mat_entry(S, i, i), expected))
        {
            flint_printf("FAIL (hardcoded):\n");
            flint_printf("S[%wd,%wd] = %{fmpz}, expected %s\n",
                i, i, fmpz_mat_entry(S, i, i), expected_snf[i]);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check S is in SNF */
    if (!fmpz_mat_is_in_snf(S))
    {
        flint_printf("FAIL (hardcoded): S not in SNF\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    /* Check U * A * V == S */
    fmpz_mat_mul(T1, U, A);
    fmpz_mat_mul(T2, T1, V);
    if (!fmpz_mat_equal(T2, S))
    {
        flint_printf("FAIL (hardcoded): U*A*V != S\n");
        fflush(stdout);
        flint_abort();
    }

    /* Check |det(U)| == 1 */
    if (m > 0)
    {
        fmpz_mat_det(det, U);
        fmpz_abs(det, det);
        if (!fmpz_is_one(det))
        {
            flint_printf("FAIL (hardcoded): |det(U)| != 1\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check |det(V)| == 1 */
    if (n > 0)
    {
        fmpz_mat_det(det, V);
        fmpz_abs(det, det);
        if (!fmpz_is_one(det))
        {
            flint_printf("FAIL (hardcoded): |det(V)| != 1\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(expected);
    fmpz_clear(det);
    fmpz_mat_clear(T2);
    fmpz_mat_clear(T1);
    fmpz_mat_clear(V);
    fmpz_mat_clear(U);
    fmpz_mat_clear(S);
    fmpz_mat_clear(A);
}

TEST_FUNCTION_START(fmpz_mat_snf_transform, state)
{
    slong iter;

    /* Magma-verified hardcoded test cases */

    /* structured 5x5 */
    {
        slong data[] = {2,4,6,8,10,1,3,5,7,9,0,2,4,6,8,1,1,1,1,1,3,6,9,12,15};
        const char * snf[] = {"1","1","0","0","0"};
        _test_transform_case(5, 5, data, snf);
    }
    /* structured tall 6x4 */
    {
        slong data[] = {6,12,18,24,3,9,15,21,0,6,12,18,3,3,3,3,9,18,27,36,12,24,36,48};
        const char * snf[] = {"3","3","0","0"};
        _test_transform_case(6, 4, data, snf);
    }
    /* structured wide 3x7 */
    {
        slong data[] = {10,20,30,40,50,60,70,5,15,25,35,45,55,65,0,10,20,30,40,50,60};
        const char * snf[] = {"5","10","0"};
        _test_transform_case(3, 7, data, snf);
    }
    /* random square 10x10 */
    {
        slong data[] = {5,11,-18,15,-14,7,-15,5,-11,11,-12,-8,7,9,16,-4,-10,16,-12,17,17,0,-8,-4,20,16,-5,6,-15,8,6,-2,-1,-6,2,20,8,-8,12,-19,11,-14,7,20,-16,14,-20,11,-7,-6,18,-17,9,6,-5,8,16,9,11,3,-3,-12,10,0,6,14,17,4,13,5,-2,-16,-10,6,-13,-9,13,19,1,-5,-11,-11,3,-16,17,6,-15,-7,18,16,19,-12,-5,-17,-15,19,-5,2,2,15};
        const char * snf[] = {"1","1","1","1","1","1","1","2","4","3395595701880"};
        _test_transform_case(10, 10, data, snf);
    }
    /* prescribed divisors 4x6 */
    {
        slong data[] = {0,-4466,566,6304,-2618,104,0,12894,-1634,-18200,7558,-300,0,-466,78,776,-458,76,0,11970,-1590,-17360,7754,-536};
        const char * snf[] = {"2","4","16","64"};
        _test_transform_case(4, 6, data, snf);
    }
    /* prescribed divisors 5x5 */
    {
        slong data[] = {-129276,-19962,407454,278226,31434,-2556756,-391248,8102508,5524632,632040,5358948,818754,-16998894,-11587626,-1328562,1529784,233520,-4855032,-3309072,-379872,-3692352,-562932,11726772,7991124,918996};
        const char * snf[] = {"6","12","48","144","720"};
        _test_transform_case(5, 5, data, snf);
    }
    /* prescribed divisors 4x7 */
    {
        slong data[] = {1150,1580,-7800,69820,7800,31800,2400,-2645,-3635,17940,-160585,-17940,-73140,-5520,1375,2705,-9360,81295,9360,37440,2160,725,1485,-4940,42755,4940,19720,1100};
        const char * snf[] = {"5","10","20","60"};
        _test_transform_case(4, 7, data, snf);
    }
    /* single row 1x6 */
    {
        slong data[] = {12,-18,30,0,6,-24};
        const char * snf[] = {"6"};
        _test_transform_case(1, 6, data, snf);
    }
    /* single column 7x1 */
    {
        slong data[] = {20,-30,10,0,40,-50,60};
        const char * snf[] = {"10"};
        _test_transform_case(7, 1, data, snf);
    }
    /* negative diagonal 3x3 */
    {
        slong data[] = {-6,0,0,0,-12,0,0,0,3};
        const char * snf[] = {"3","6","12"};
        _test_transform_case(3, 3, data, snf);
    }
    /* singular 8x8 rank 6 */
    {
        slong data[] = {156,-77,-31,58,144,35,-45,-23,196,-222,9,116,55,69,200,-77,-22,-9,-229,-53,-83,115,-84,178,-63,66,-98,16,-20,-62,-124,2,-94,-43,-31,-51,-59,17,45,62,-30,-23,-91,5,13,3,-79,59,18,92,191,-48,53,-25,25,-67,122,-63,84,14,-45,96,194,-25};
        const char * snf[] = {"1","1","1","1","1","1","0","0"};
        _test_transform_case(8, 8, data, snf);
    }
    /* singular 7x7 rank 1 */
    {
        slong data[] = {-30,-27,-24,-18,-6,0,-9,-60,-54,-48,-36,-12,0,-18,-10,-9,-8,-6,-2,0,-3,20,18,16,12,4,0,6,20,18,16,12,4,0,6,-30,-27,-24,-18,-6,0,-9,-30,-27,-24,-18,-6,0,-9};
        const char * snf[] = {"1","0","0","0","0","0","0"};
        _test_transform_case(7, 7, data, snf);
    }

    /* Randomized tests: square, non-square, rank-deficient */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S2, U, V, T1, T2;
        fmpz_t det;
        slong m, n, r, b, d, i;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        r = n_randint(state, FLINT_MIN(m, n) + 1);
        d = FLINT_MIN(m, n);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S2, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(V, n, n);
        fmpz_mat_init(T1, m, n);
        fmpz_mat_init(T2, m, n);
        fmpz_init(det);

        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, r, b);

        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, n_randint(state, 2 * m * n + 1));

        fmpz_mat_snf_transform(S, U, V, A);

        /* Check S is in SNF */
        if (!fmpz_mat_is_in_snf(S))
        {
            flint_printf("FAIL:\n");
            flint_printf("S not in SNF, m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Check U * A * V == S */
        fmpz_mat_mul(T1, U, A);
        fmpz_mat_mul(T2, T1, V);
        if (!fmpz_mat_equal(T2, S))
        {
            flint_printf("FAIL:\n");
            flint_printf("U*A*V != S, m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Check |det(U)| == 1 (unimodular) */
        if (m > 0)
        {
            fmpz_mat_det(det, U);
            fmpz_abs(det, det);
            if (!fmpz_is_one(det))
            {
                flint_printf("FAIL:\n");
                flint_printf("|det(U)| != 1, m=%wd n=%wd\n", m, n);
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check |det(V)| == 1 (unimodular) */
        if (n > 0)
        {
            fmpz_mat_det(det, V);
            fmpz_abs(det, det);
            if (!fmpz_is_one(det))
            {
                flint_printf("FAIL:\n");
                flint_printf("|det(V)| != 1, m=%wd n=%wd\n", m, n);
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check rank: number of nonzero diagonal entries */
        {
            slong snf_rank = 0;
            for (i = 0; i < d; i++)
                if (!fmpz_is_zero(fmpz_mat_entry(S, i, i)))
                    snf_rank++;

            if (snf_rank != fmpz_mat_rank(A))
            {
                flint_printf("FAIL:\n");
                flint_printf("snf rank %wd != matrix rank %wd\n",
                    snf_rank, fmpz_mat_rank(A));
                fflush(stdout);
                flint_abort();
            }
        }

        /* Cross-check against fmpz_mat_snf */
        fmpz_mat_snf(S2, A);
        if (!fmpz_mat_equal(S, S2))
        {
            flint_printf("FAIL:\n");
            flint_printf("snf_transform disagrees with snf\n");
            flint_printf("m=%wd n=%wd r=%wd\n", m, n, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            flint_printf("snf_transform: ");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            flint_printf("snf:           ");
            fmpz_mat_print_pretty(S2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(det);
        fmpz_mat_clear(T2);
        fmpz_mat_clear(T1);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S2);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    /* Edge cases */
    {
        /* 0x0 */
        fmpz_mat_t S, U, V;
        fmpz_mat_init(S, 0, 0);
        fmpz_mat_init(U, 0, 0);
        fmpz_mat_init(V, 0, 0);
        fmpz_mat_snf_transform(S, U, V, S);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
    }
    {
        /* 0x5 */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 0, 5);
        fmpz_mat_init(S, 0, 5);
        fmpz_mat_init(U, 0, 0);
        fmpz_mat_init(V, 5, 5);
        fmpz_mat_snf_transform(S, U, V, A);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* 3x0 */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 3, 0);
        fmpz_mat_init(S, 3, 0);
        fmpz_mat_init(U, 3, 3);
        fmpz_mat_init(V, 0, 0);
        fmpz_mat_snf_transform(S, U, V, A);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* 1x1 negative */
        fmpz_mat_t A, S, U, V, T1, T2;
        fmpz_mat_init(A, 1, 1);
        fmpz_mat_init(S, 1, 1);
        fmpz_mat_init(U, 1, 1);
        fmpz_mat_init(V, 1, 1);
        fmpz_mat_init(T1, 1, 1);
        fmpz_mat_init(T2, 1, 1);

        fmpz_set_si(fmpz_mat_entry(A, 0, 0), -7);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_equal_si(fmpz_mat_entry(S, 0, 0), 7))
        {
            flint_printf("FAIL: 1x1 negative\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_mul(T1, U, A);
        fmpz_mat_mul(T2, T1, V);
        if (!fmpz_mat_equal(T2, S))
        {
            flint_printf("FAIL: U*A*V != S for 1x1\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(T2);
        fmpz_mat_clear(T1);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* Identity matrix */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 5, 5);
        fmpz_mat_init(S, 5, 5);
        fmpz_mat_init(U, 5, 5);
        fmpz_mat_init(V, 5, 5);

        fmpz_mat_one(A);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_mat_is_one(S))
        {
            flint_printf("FAIL: identity\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* Zero matrix */
        fmpz_mat_t A, S, U, V;
        fmpz_mat_init(A, 4, 6);
        fmpz_mat_init(S, 4, 6);
        fmpz_mat_init(U, 4, 4);
        fmpz_mat_init(V, 6, 6);

        fmpz_mat_zero(A);
        fmpz_mat_snf_transform(S, U, V, A);

        if (!fmpz_mat_is_zero(S))
        {
            flint_printf("FAIL: zero matrix\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    /* NULL transform tests */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S, S_full, U, U_full, V, V_full;
        slong m, n, r, b;

        m = n_randint(state, 15);
        n = n_randint(state, 15);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);
        fmpz_mat_init(S_full, m, n);
        fmpz_mat_init(U, m, m);
        fmpz_mat_init(U_full, m, m);
        fmpz_mat_init(V, n, n);
        fmpz_mat_init(V_full, n, n);
        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, r, b);
        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, n_randint(state, 2 * m * n + 1));

        /* Full transform for reference */
        fmpz_mat_snf_transform(S_full, U_full, V_full, A);

        /* NULL V: only U computed */
        fmpz_mat_snf_transform(S, U, NULL, A);
        if (!fmpz_mat_equal(S, S_full))
        {
            flint_printf("FAIL (NULL V): S differs\n");
            fflush(stdout);
            flint_abort();
        }
        if (!fmpz_mat_equal(U, U_full))
        {
            flint_printf("FAIL (NULL V): U differs\n");
            fflush(stdout);
            flint_abort();
        }

        /* NULL U: only V computed */
        fmpz_mat_snf_transform(S, NULL, V, A);
        if (!fmpz_mat_equal(S, S_full))
        {
            flint_printf("FAIL (NULL U): S differs\n");
            fflush(stdout);
            flint_abort();
        }
        if (!fmpz_mat_equal(V, V_full))
        {
            flint_printf("FAIL (NULL U): V differs\n");
            fflush(stdout);
            flint_abort();
        }

        /* Both NULL: just SNF */
        fmpz_mat_snf_transform(S, NULL, NULL, A);
        if (!fmpz_mat_equal(S, S_full))
        {
            flint_printf("FAIL (both NULL): S differs\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(V_full);
        fmpz_mat_clear(U_full);
        fmpz_mat_clear(V);
        fmpz_mat_clear(U);
        fmpz_mat_clear(S_full);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}

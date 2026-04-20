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
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/*
    Helper: build m x n matrix from slong array, compute elementary_divisors,
    and check against expected values (verified with Magma).
*/
static void
_test_ed_case(slong m, slong n, slong expected_rank,
    const slong * data, const slong * expected_ed)
{
    fmpz_mat_t A;
    fmpz * ed;
    slong rank, d, i, j;

    d = FLINT_MIN(m, n);

    fmpz_mat_init(A, m, n);
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), data[i * n + j]);

    ed = _fmpz_vec_init(d);
    rank = fmpz_mat_elementary_divisors(ed, A);

    if (rank != expected_rank)
    {
        flint_printf("FAIL (hardcoded): rank %wd != expected %wd\n",
            rank, expected_rank);
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < expected_rank; i++)
    {
        if (!fmpz_equal_si(&ed[i], expected_ed[i]))
        {
            flint_printf("FAIL (hardcoded): ed[%wd] = %{fmpz}, "
                "expected %wd\n", i, &ed[i], expected_ed[i]);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }
    }

    _fmpz_vec_clear(ed, d);
    fmpz_mat_clear(A);
}

TEST_FUNCTION_START(fmpz_mat_elementary_divisors, state)
{
    slong iter;

    /* Magma-verified hardcoded test cases */

    /* structured 5x5 rank 2 */
    {
        slong data[] = {2,4,6,8,10,1,3,5,7,9,0,2,4,6,8,1,1,1,1,1,3,6,9,12,15};
        slong ed[] = {1,1};
        _test_ed_case(5, 5, 2, data, ed);
    }
    /* structured tall 6x4 rank 2 */
    {
        slong data[] = {6,12,18,24,3,9,15,21,0,6,12,18,3,3,3,3,9,18,27,36,12,24,36,48};
        slong ed[] = {3,3};
        _test_ed_case(6, 4, 2, data, ed);
    }
    /* structured wide 3x7 rank 2 */
    {
        slong data[] = {10,20,30,40,50,60,70,5,15,25,35,45,55,65,0,10,20,30,40,50,60};
        slong ed[] = {5,10};
        _test_ed_case(3, 7, 2, data, ed);
    }
    /* random square 7x7 rank 7 */
    {
        slong data[] = {-19,2,-1,-4,-20,15,-18,16,-9,-15,11,-17,-2,14,8,8,6,14,-10,-5,11,-16,20,-2,-2,-18,-10,1,20,-20,6,8,-17,-3,-18,7,-7,-4,-12,9,-2,-18,-1,-10,-4,-2,-14,-16,-1};
        slong ed[] = {1,1,1,1,1,1,3302157};
        _test_ed_case(7, 7, 7, data, ed);
    }
    /* random square 5x5 rank 5 */
    {
        slong data[] = {3,8,-10,-12,-15,7,17,0,10,1,3,18,-13,10,-5,4,7,-8,15,-17,-10,6,-1,-3,1};
        slong ed[] = {1,1,1,1,899640};
        _test_ed_case(5, 5, 5, data, ed);
    }
    /* prescribed 7x7 */
    {
        slong data[] = {25699,29914,1055,-4872,2376,16678,7213,1677312,2047248,77120,-126256,64080,1329120,572880,-1026034,-1247384,-46702,87828,-44268,-799656,-344746,1621804,2005872,76780,-67632,36112,1353576,582916,420416,513044,19316,-31916,16196,332804,143444,-320188,-402056,-15652,1080,-1320,-282672,-121612,-1274215,-1512490,-54819,185272,-91224,-897246,-387481};
        slong ed[] = {1,4,8,40,160,480,1920};
        _test_ed_case(7, 7, 7, data, ed);
    }
    /* prescribed 7x7 (b) */
    {
        slong data[] = {68026,-52588,-16032,4004,7008,-7506,17010,4475657,-2220856,-514916,127678,223356,-417867,926745,221217,-180656,-48184,11918,20836,-30027,53745,-825828,349624,80500,-20012,-35024,66168,-165980,391273,-168544,-56184,14042,24584,-21063,85105,-50000,52540,16016,-4000,-7000,7500,-14000,583694,-325572,-72440,17876,31252,-63114,122590};
        slong ed[] = {1,4,20,100,500,1500,3000};
        _test_ed_case(7, 7, 7, data, ed);
    }
    /* singular 12x12 rank 6 */
    {
        slong data[] = {9,76,-49,-98,50,138,-95,39,-132,89,21,70,25,155,22,-36,-45,-110,-61,23,-63,-116,45,99,210,112,-29,-21,-120,-57,-36,13,24,-96,88,74,-117,24,67,223,-21,-51,87,-58,51,-43,-114,-2,-53,66,14,-93,34,13,-77,26,-106,31,23,38,191,-22,67,-60,-187,-163,-14,-31,105,-41,121,-94,-167,-3,-59,25,16,70,-84,-8,-163,14,-132,-37,-14,-47,6,-53,38,-62,25,17,38,-41,30,-13,-29,104,-31,-172,23,-11,-140,55,-170,-27,46,66,123,83,85,-28,-98,-156,5,-5,73,-74,111,31,3,10,-74,14,13,0,-5,19,-22,-87,-33,36,13,-39,3,171,-66,34,64,-56,70,18,-79,-51};
        slong ed[] = {1,1,1,1,1,1};
        _test_ed_case(12, 12, 6, data, ed);
    }

    /* Randomized tests: cross-check against SNF diagonal */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, S;
        fmpz * ed;
        slong m, n, r, b, d, i, rank;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        r = n_randint(state, FLINT_MIN(m, n) + 1);
        d = FLINT_MIN(m, n);

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(S, m, n);

        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        fmpz_mat_randrank(A, state, r, b);

        if (n_randint(state, 2))
            fmpz_mat_randops(A, state,
                n_randint(state, 2 * m * n + 1));

        ed = _fmpz_vec_init(d);
        rank = fmpz_mat_elementary_divisors(ed, A);

        /* Check rank */
        if (rank != r)
        {
            flint_printf("FAIL:\n");
            flint_printf("rank %wd != expected %wd\n", rank, r);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        /* Cross-check against SNF diagonal */
        fmpz_mat_snf(S, A);

        for (i = 0; i < r; i++)
        {
            if (!fmpz_equal(&ed[i], fmpz_mat_entry(S, i, i)))
            {
                flint_printf("FAIL:\n");
                flint_printf("ed[%wd] = %{fmpz} != snf[%wd] = %{fmpz}\n",
                    i, &ed[i], i, fmpz_mat_entry(S, i, i));
                flint_printf("m=%wd n=%wd r=%wd\n", m, n, r);
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpz_vec_clear(ed, d);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }

    /* Edge cases */
    {
        /* 0x0 matrix */
        fmpz_mat_t A;
        fmpz_mat_init(A, 0, 0);
        if (fmpz_mat_elementary_divisors(NULL, A) != 0)
        {
            flint_printf("FAIL: 0x0 rank\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }
    {
        /* Zero matrix */
        fmpz_mat_t A;
        fmpz_mat_init(A, 5, 3);
        fmpz_mat_zero(A);
        if (fmpz_mat_elementary_divisors(NULL, A) != 0)
        {
            flint_printf("FAIL: zero matrix rank\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mat_clear(A);
    }
    {
        /* Identity matrix */
        fmpz_mat_t A;
        fmpz * ed;
        slong i;
        fmpz_mat_init(A, 4, 4);
        fmpz_mat_one(A);
        ed = _fmpz_vec_init(4);
        if (fmpz_mat_elementary_divisors(ed, A) != 4)
        {
            flint_printf("FAIL: identity rank\n");
            fflush(stdout);
            flint_abort();
        }
        for (i = 0; i < 4; i++)
        {
            if (!fmpz_is_one(&ed[i]))
            {
                flint_printf("FAIL: identity ed[%wd] = %{fmpz}\n",
                    i, &ed[i]);
                fflush(stdout);
                flint_abort();
            }
        }
        _fmpz_vec_clear(ed, 4);
        fmpz_mat_clear(A);
    }
    {
        /* Known diagonal matrix */
        fmpz_mat_t A;
        fmpz * ed;
        fmpz_mat_init(A, 3, 3);
        fmpz_mat_zero(A);
        fmpz_set_si(fmpz_mat_entry(A, 0, 0), 12);
        fmpz_set_si(fmpz_mat_entry(A, 1, 1), 6);
        fmpz_set_si(fmpz_mat_entry(A, 2, 2), 0);
        ed = _fmpz_vec_init(3);
        if (fmpz_mat_elementary_divisors(ed, A) != 2)
        {
            flint_printf("FAIL: diag rank\n");
            fflush(stdout);
            flint_abort();
        }
        if (!fmpz_equal_si(&ed[0], 6)
            || !fmpz_equal_si(&ed[1], 12))
        {
            flint_printf("FAIL: diag ed = [%{fmpz}, %{fmpz}], "
                "expected [6, 12]\n", &ed[0], &ed[1]);
            fflush(stdout);
            flint_abort();
        }
        _fmpz_vec_clear(ed, 3);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}

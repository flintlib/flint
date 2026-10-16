/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "ulong_extras.h"

/*
    Exercises code paths of the floating-point LLL variants not covered by
    the other tests: inputs with entries too large to use the packed
    representation, inputs that make the double precision Babai step fail
    (a huge entry next to tiny ones underflows in the scaled double
    approximation), and linearly dependent rows in the multiprecision LLL,
    with packed and unpacked matrices.

    Since these inputs are ill-conditioned, the output is not required to
    be reduced; what is verified is that the output is the given unimodular
    transformation of the input, which checks that all row operations are
    performed exactly.
*/

/* random small matrix with one entry of about hugebits bits */
static void
randmat_huge(fmpz_mat_t A, flint_rand_t state, slong hugebits)
{
    fmpz * c;
    fmpz_mat_randtest(A, state, 1 + n_randint(state, 30));
    c = fmpz_mat_entry(A, n_randint(state, A->r), n_randint(state, A->c));
    fmpz_randbits(c, state, hugebits);
    if (n_randint(state, 2))
        fmpz_neg(c, c);
}

/* is B = U * A with U unimodular? */
static int
check_transform(const fmpz_mat_t A, const fmpz_mat_t U, const fmpz_mat_t B)
{
    fmpz_mat_t T;
    fmpz_t det;
    int ok;

    fmpz_mat_init(T, B->r, B->c);
    fmpz_init(det);
    fmpz_mat_mul(T, U, A);
    fmpz_mat_det(det, U);
    ok = fmpz_mat_equal(T, B) && fmpz_is_pm1(det);
    fmpz_mat_clear(T);
    fmpz_clear(det);
    return ok;
}

TEST_FUNCTION_START(fmpz_lll_unpacked, state)
{
    slong iter;

    /* double precision LLL: unpacked (huge entries) and packed but failing (underflow) */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, B, U;
        fmpz_lll_t fl;
        fmpz_t gs_B;
        slong d, n, hugebits;
        const slong packed_max_bits = FMPZ_LLL_PACKED_MAX_LIMBS * FLINT_BITS;
        int res, kind = n_randint(state, 4), unpacked = n_randint(state, 2);

        d = 2 + n_randint(state, 8);
        n = d + n_randint(state, 3);

        /*
            Packed, but as large as possible: with 64-bit limbs the entry
            exceeds the double exponent range relative to the small entries,
            which makes the Babai step fail (with 32-bit limbs the packable
            range is too small for this to happen).
        */
        if (unpacked)
            hugebits = FMPZ_LLL_PACKED_MAX_LIMBS * FLINT_BITS + 100 + n_randint(state, 200);
        else
            hugebits = packed_max_bits * 5 / 8 + n_randint(state, packed_max_bits * 3 / 8 - 64);

        fmpz_mat_init(A, d, n);
        fmpz_mat_init(B, d, n);
        fmpz_mat_init(U, d, d);
        fmpz_mat_one(U);
        fmpz_init(gs_B);
        randmat_huge(A, state, hugebits);
        fmpz_mat_set(B, A);
        fmpz_set_ui(gs_B, 1);
        fmpz_mul_2exp(gs_B, gs_B, n_randint(state, 200));
        fmpz_lll_context_init(fl, 0.75, 0.81, Z_BASIS, APPROX);

        if (unpacked != (fmpz_lll_packed_limbs(A) > FMPZ_LLL_PACKED_MAX_LIMBS))
            TEST_FUNCTION_FAIL("unexpected packing decision\n");

        if (kind == 0)
            res = fmpz_lll_d(B, U, fl);
        else if (kind == 1)
            res = fmpz_lll_d_heuristic(B, U, fl);
        else if (kind == 2)
            res = fmpz_lll_d_with_removal(B, U, gs_B, fl);
        else
            res = fmpz_lll_d_with_removal_knapsack(B, U, gs_B, fl);

        if (!check_transform(A, U, B))
            TEST_FUNCTION_FAIL("lll_d kind %d res %d: output is not a unimodular transformation of the input\nA = %{fmpz_mat}\nB = %{fmpz_mat}\nU = %{fmpz_mat}\n", kind, res, A, B, U);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(U);
        fmpz_clear(gs_B);
    }

    /* knapsack variant with early size reduction of the later rows, with Babai failures */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, B, U;
        fmpz_lll_t fl;
        fmpz_t gs_B;
        slong d, n, i;
        const slong packed_max_bits = FMPZ_LLL_PACKED_MAX_LIMBS * FLINT_BITS;
        int res;

        d = 6 + n_randint(state, 8);
        n = d + n_randint(state, 3);

        fmpz_mat_init(A, d, n);
        fmpz_mat_init(B, d, n);
        fmpz_mat_init(U, d, d);
        fmpz_mat_one(U);
        fmpz_init(gs_B);
        fmpz_mat_randtest(A, state, 1 + n_randint(state, 30));

        /* several rows with a huge entry each */
        for (i = 0; i < d; i++)
            if (n_randint(state, 3) == 0)
                fmpz_randbits(fmpz_mat_entry(A, i, n_randint(state, n)), state,
                    packed_max_bits * 5 / 8 + n_randint(state, packed_max_bits * 3 / 8 - 64));

        fmpz_mat_set(B, A);
        fmpz_set_ui(gs_B, 1);
        fmpz_mul_2exp(gs_B, gs_B, n_randint(state, 200));
        fmpz_lll_context_init(fl, 0.75, 0.81, Z_BASIS, APPROX);

        res = fmpz_lll_d_with_removal_knapsack(B, U, gs_B, fl);

        if (!check_transform(A, U, B))
            TEST_FUNCTION_FAIL("knapsack res %d: output is not a unimodular transformation of the input\nA = %{fmpz_mat}\nB = %{fmpz_mat}\nU = %{fmpz_mat}\n", res, A, B, U);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(U);
        fmpz_clear(gs_B);
    }

    /* multiprecision LLL: unpacked (m > FMPZ_LLL_PACKED_MAX_LIMBS_MPF) and linearly dependent rows */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A, B, U;
        fmpz_lll_t fl;
        fmpz_t gs_B;
        slong d, n, i;
        int res, with_removal = n_randint(state, 2), unpacked = n_randint(state, 2);
        flint_bitcnt_t prec = 64 * (1 + n_randint(state, 3));

        d = 2 + n_randint(state, 8);
        n = 1 + n_randint(state, d);   /* usually more rows than columns */

        fmpz_mat_init(A, d, n);
        fmpz_mat_init(B, d, n);
        fmpz_mat_init(U, d, d);
        fmpz_mat_one(U);
        fmpz_init(gs_B);

        fmpz_mat_randtest(A, state, 1 + n_randint(state, 60));

        /* duplicated row */
        if (n_randint(state, 2))
        {
            i = n_randint(state, d - 1);
            _fmpz_vec_set(fmpz_mat_row(A, i + 1), fmpz_mat_row(A, i), n);
        }

        if (unpacked)
        {
            fmpz * c = fmpz_mat_entry(A, n_randint(state, d), n_randint(state, n));
            fmpz_randbits(c, state, FMPZ_LLL_PACKED_MAX_LIMBS_MPF * FLINT_BITS + 50 + n_randint(state, 300));
        }

        fmpz_mat_set(B, A);
        fmpz_set_ui(gs_B, 1);
        fmpz_mul_2exp(gs_B, gs_B, n_randint(state, 50));
        fmpz_lll_context_init(fl, 0.75, 0.81, Z_BASIS, APPROX);

        if (unpacked != (fmpz_lll_packed_limbs(A) > FMPZ_LLL_PACKED_MAX_LIMBS_MPF))
            TEST_FUNCTION_FAIL("unexpected packing decision\n");

        if (with_removal)
            res = fmpz_lll_mpf2_with_removal(B, U, prec, gs_B, fl);
        else
            res = fmpz_lll_mpf2(B, U, prec, fl);

        if (!check_transform(A, U, B))
            TEST_FUNCTION_FAIL("mpf2 removal %d res %d: output is not a unimodular transformation of the input\nA = %{fmpz_mat}\nB = %{fmpz_mat}\nU = %{fmpz_mat}\n", with_removal, res, A, B, U);

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(U);
        fmpz_clear(gs_B);
    }

    TEST_FUNCTION_END(state);
}

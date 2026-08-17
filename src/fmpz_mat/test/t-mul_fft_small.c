/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fft_small.h"
#include "fmpz_mat.h"

/* Transformed-integer matrix product against the default algorithm. The
   function returns 0 when the representation is unavailable or judged
   unprofitable; entries are drawn dense (fmpz_mat_randtest's sparse
   entries make degenerate shapes) with both signs; every fourth
   iteration squares (B aliased to A), and a rotating storage budget
   forces each blocking tier of the driver. */
TEST_FUNCTION_START(fmpz_mat_mul_fft_small, state)
{
    for (slong iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong r = 1 + n_randint(state, 4);
        slong k = 1 + n_randint(state, 4);
        slong c = 1 + n_randint(state, 4);

        /* asymmetric shapes reach the one-side-resident streaming tiers
           in both orientations, and a large inner dimension with a
           starved budget reaches the inner-blocked accumulation tier */
        if (iter % 8 == 5)
        { r = 1 + n_randint(state, 2); c = 5 + n_randint(state, 4); }
        else if (iter % 8 == 6)
        { r = 5 + n_randint(state, 4); c = 1 + n_randint(state, 2); }
        else if (iter % 8 == 7)
        { r = 1; c = 1; k = 8 + n_randint(state, 5); }
        ulong save_limit = flint_fft_small_max_transformed_ring_size;
        int square = (iter % 4 == 1);

        if (square)
        { k = r; c = r; }
        /* rotate the budget to force the resident, streamed, doubly
           blocked and inner-blocked strategies in turn */
        if (iter % 4 == 2)
            flint_fft_small_max_transformed_ring_size = UWORD(1) << 20;
        else if (iter % 4 == 3)
            flint_fft_small_max_transformed_ring_size =
                UWORD(1) << (13 + n_randint(state, 8));
        slong bits = 8192 + n_randint(state, 8192);
        fmpz_mat_t A, B, C, D, A0, B0;
        slong i, j;

        fmpz_mat_init(A, r, k);
        fmpz_mat_init(B, k, c);
        fmpz_mat_init(C, r, c);
        fmpz_mat_init(D, r, c);

        for (i = 0; i < r; i++)
            for (j = 0; j < k; j++)
            {
                fmpz_randbits(fmpz_mat_entry(A, i, j), state, bits);
                /* every third iteration keeps everything nonnegative to
                   exercise the unsigned-context fast path, which signed
                   entries almost never select by chance */
                if (iter % 3 == 0)
                    fmpz_abs(fmpz_mat_entry(A, i, j),
                             fmpz_mat_entry(A, i, j));
                if (n_randint(state, 8) == 0)
                    fmpz_zero(fmpz_mat_entry(A, i, j));
            }
        for (i = 0; i < k; i++)
            for (j = 0; j < c; j++)
            {
                fmpz_randbits(fmpz_mat_entry(B, i, j), state, bits);
                if (iter % 3 == 0)
                    fmpz_abs(fmpz_mat_entry(B, i, j),
                             fmpz_mat_entry(B, i, j));
                if (n_randint(state, 8) == 0)
                    fmpz_zero(fmpz_mat_entry(B, i, j));
            }

        fmpz_mat_init_set(A0, A);
        fmpz_mat_init_set(B0, B);

        fmpz_mat_mul(D, A, square ? A : B);

        if (fmpz_mat_mul_fft_small(C, A, square ? A : B))
        {
            if (!fmpz_mat_equal(C, D))
                TEST_FUNCTION_FAIL(
                    "wrong product\n"
                    "dims %wd x %wd x %wd, entry bits %wd\n",
                    r, k, c, bits);
        }

        /* the inputs must be untouched whether or not the product ran */
        if (!fmpz_mat_equal(A, A0) || !fmpz_mat_equal(B, B0))
            TEST_FUNCTION_FAIL("input matrices modified\n");

        /* truncated product: each entry's magnitude above limb lo, with
           the entry's sign and at most one unit of error in the lowest
           returned limb */
        {
            slong lo = n_randint(state, 2 + bits / (2 * FLINT_BITS));
            fmpz_mat_t E;

            fmpz_mat_init(E, r, c);
            if (fmpz_mat_mul_fft_small_trunc(E, A, square ? A : B, lo))
            {
                for (i = 0; i < r; i++)
                    for (j = 0; j < c; j++)
                    {
                        fmpz_t want, got, diff;

                        fmpz_init(want);
                        fmpz_init(got);
                        fmpz_init(diff);
                        fmpz_abs(want, fmpz_mat_entry(D, i, j));
                        fmpz_fdiv_q_2exp(want, want, FLINT_BITS * lo);
                        fmpz_abs(got, fmpz_mat_entry(E, i, j));
                        fmpz_sub(diff, got, want);
                        fmpz_abs(diff, diff);

                        if (fmpz_cmp_ui(diff, 1) > 0)
                            TEST_FUNCTION_FAIL(
                                "truncated entry (%wd, %wd) off by more "
                                "than one unit\nlo = %wd, bits = %wd\n",
                                i, j, lo, bits);

                        if (!fmpz_is_zero(got) && !fmpz_is_zero(want) &&
                                fmpz_sgn(fmpz_mat_entry(E, i, j)) !=
                                fmpz_sgn(fmpz_mat_entry(D, i, j)))
                            TEST_FUNCTION_FAIL(
                                "truncated entry (%wd, %wd) sign wrong\n",
                                i, j);

                        fmpz_clear(want);
                        fmpz_clear(got);
                        fmpz_clear(diff);
                    }
            }
            fmpz_mat_clear(E);
        }

        flint_fft_small_max_transformed_ring_size = save_limit;
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(A0);
        fmpz_mat_clear(B0);
    }

        /* structured matrices: entries from {0, +1, -1, big} */
    for (slong it2 = 0; it2 < 10 * flint_test_multiplier(); it2++)
    {
        slong r = 1 + n_randint(state, 6), c = 1 + n_randint(state, 6),
              k = 1 + n_randint(state, 6);
        slong bits = 11000 + n_randint(state, 6000);
        fmpz_mat_t A, B, C, D;
        fmpz_mat_init(A, r, k); fmpz_mat_init(B, k, c);
        fmpz_mat_init(C, r, c); fmpz_mat_init(D, r, c);
        for (slong i2 = 0; i2 < r; i2++)
            for (slong j2 = 0; j2 < k; j2++)
                switch (n_randint(state, 4))
                {
                    case 0: break;
                    case 1: fmpz_set_si(fmpz_mat_entry(A, i2, j2),
                                n_randint(state, 2) ? 1 : -1); break;
                    default: fmpz_randtest(fmpz_mat_entry(A, i2, j2),
                                           state, bits);
                }
        for (slong i2 = 0; i2 < k; i2++)
            for (slong j2 = 0; j2 < c; j2++)
                switch (n_randint(state, 4))
                {
                    case 0: break;
                    case 1: fmpz_set_si(fmpz_mat_entry(B, i2, j2),
                                n_randint(state, 2) ? 1 : -1); break;
                    default: fmpz_randtest(fmpz_mat_entry(B, i2, j2),
                                           state, bits);
                }
        if (fmpz_mat_mul_fft_small(C, A, B))
        {
            fmpz_mat_mul_classical(D, A, B);
            if (!fmpz_mat_equal(C, D))
                TEST_FUNCTION_FAIL("structured matrix mismatch: "
                    "%wd x %wd x %wd, bits = %wd\n", r, k, c, bits);
        }
        fmpz_mat_clear(A); fmpz_mat_clear(B);
        fmpz_mat_clear(C); fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}

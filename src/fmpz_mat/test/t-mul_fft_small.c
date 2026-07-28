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
#include "fmpz_mat.h"

/* Transformed-integer matrix product against the default algorithm. The
   function returns 0 when the representation is unavailable or judged
   unprofitable; entries are drawn dense (fmpz_mat_randtest is sparse
   enough to fall below the gate) with both signs. */
TEST_FUNCTION_START(fmpz_mat_mul_fft_small, state)
{
    for (slong iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong r = 1 + n_randint(state, 4);
        slong k = 1 + n_randint(state, 4);
        slong c = 1 + n_randint(state, 4);
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
                if (n_randint(state, 8) == 0)
                    fmpz_zero(fmpz_mat_entry(A, i, j));
            }
        for (i = 0; i < k; i++)
            for (j = 0; j < c; j++)
            {
                fmpz_randbits(fmpz_mat_entry(B, i, j), state, bits);
                if (n_randint(state, 8) == 0)
                    fmpz_zero(fmpz_mat_entry(B, i, j));
            }

        fmpz_mat_init_set(A0, A);
        fmpz_mat_init_set(B0, B);

        fmpz_mat_mul(D, A, B);

        if (fmpz_mat_mul_fft_small(C, A, B))
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
            if (fmpz_mat_mul_fft_small_trunc(E, A, B, lo))
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

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
        fmpz_mat_clear(A0);
        fmpz_mat_clear(B0);
    }

    TEST_FUNCTION_END(state);
}

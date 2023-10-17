/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul, state)
{
    fmpz_mat_t A, B, C, D;
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong m, n, k;
        slong abits, bbits, bits;

        if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 50);
            n = n_randint(state, 50);
            k = n_randint(state, 50);
        }
        else
        {
            m = n_randint(state, 8);
            n = n_randint(state, 8);
            k = n_randint(state, 8);
        }

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        abits = fmpz_mat_max_bits(A);
        bbits = fmpz_mat_max_bits(B);
        abits = FLINT_ABS(abits);
        bbits = FLINT_ABS(bbits);
        bits = abits + bbits + FLINT_BIT_COUNT(n) + 1;

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul_classical_inline(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n\n");
            fmpz_mat_print(A); flint_printf("\n\n");
            fmpz_mat_print(B); flint_printf("\n\n");
            fmpz_mat_print(C); flint_printf("\n\n");
            fmpz_mat_print(D); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (bits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_mat_mul_small_1(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_small_1)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (abits <= SMALL_FMPZ_BITCOUNT_MAX && bbits <= SMALL_FMPZ_BITCOUNT_MAX && bits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_mat_mul_small_2a(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_small_2a)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (abits <= SMALL_FMPZ_BITCOUNT_MAX && bbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_mat_mul_small_2b(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_small_2b)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (abits < 2 * FLINT_BITS && bbits < 2 * FLINT_BITS)
        {
            _fmpz_mat_mul_double_word(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_double_word)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                fflush(stdout);
                flint_abort();
            }
        }

        if (n == k)
        {
            fmpz_mat_mul(A, A, B);

            if (!fmpz_mat_equal(A, C))
            {
                flint_printf("FAIL: aliasing failed\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    /* Test aliasing with windows */
    {
        fmpz_mat_t A, B, A_window;

        fmpz_mat_init(A, 2, 2);
        fmpz_mat_init(B, 2, 2);

        fmpz_mat_window_init(A_window, A, 0, 0, 2, 2);

        fmpz_mat_one(A);
        fmpz_mat_one(B);
        fmpz_set_ui(fmpz_mat_entry(B, 0, 1), 1);
        fmpz_set_ui(fmpz_mat_entry(B, 1, 0), 1);

        fmpz_mat_mul(A_window, B, A_window);

        if (!fmpz_mat_equal(A, B))
        {
            flint_printf("FAIL: window aliasing failed\n");
	    fmpz_mat_print(A); flint_printf("\n\n");
	    fmpz_mat_print(B); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_window_clear(A_window);
        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}

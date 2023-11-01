/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"

void _fmpz_mat_full(fmpz_mat_t A, flint_bitcnt_t bits)
{
    slong i, j;

    for (i = 0; i < A->r; i++ )
    {
        for (j = 0; j < A->c; j++)
        {
            fmpz_one(fmpz_mat_entry(A, i, j));
            fmpz_mul_2exp(fmpz_mat_entry(A, i, j), fmpz_mat_entry(A, i, j), bits);
            fmpz_sub_ui(fmpz_mat_entry(A, i, j), fmpz_mat_entry(A, i, j), 1);
        }
    }
}

TEST_FUNCTION_START(fmpz_mat_mul_double_word, state)
{
    fmpz_mat_t A, B, C, D;
    slong i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int sign;
        slong m, n, k;
        flint_bitcnt_t abits, bbits;

        sign = n_randint(state, 2);
        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);
        abits = n_randint(state, 2*FLINT_BITS - sign) + 1;
        bbits = n_randint(state, 2*FLINT_BITS - sign) + 1;

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m, n);

        if (sign)
        {
            if (n_randint(state, 2))
            {
                _fmpz_mat_full(A, abits);
                _fmpz_mat_full(B, bbits);

                if (n_randint(state, 2))
                    fmpz_mat_neg(A, A);

                if (n_randint(state, 2))
                    fmpz_mat_neg(B, B);
            }
            else
            {
                fmpz_mat_randtest(A, state, abits);
                fmpz_mat_randtest(B, state, bbits);
            }
        }
        else
        {
            if (n_randint(state, 2))
            {
                _fmpz_mat_full(A, abits);
                _fmpz_mat_full(B, bbits);
            }
            else
            {
                fmpz_mat_randtest_unsigned(A, state, abits);
                fmpz_mat_randtest_unsigned(B, state, bbits);
            }
        }

        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        _fmpz_mat_mul_double_word(C, A, B);
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

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    TEST_FUNCTION_END(state);
}

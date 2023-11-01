/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(fmpz_mat_mul_fft, state)
{
    fmpz_mat_t A, B, C, D;
    slong i;

    for (i = 0; i < 100*flint_test_multiplier(); i++)
    {
        slong m, n, k;

        m = 1 + n_randint(state, 3);
        k = 1 + n_randint(state, 3);
        n = 1 + n_randint(state, 3);

        fmpz_mat_init(A, m, k);
        fmpz_mat_init(B, k, n);
        fmpz_mat_init(C, m, n);
        fmpz_mat_init(D, m, n);

        do {
            fmpz_mat_randtest(A, state, n_randint(state, 99999) + 1);
            fmpz_mat_randtest(B, state, n_randint(state, 99999) + 1);
        } while (fmpz_mat_is_zero(A) || fmpz_mat_is_zero(B));

        fmpz_mat_randtest(C, state, n_randint(state, 2000) + 1);

        fmpz_mat_mul_fft(C, A, B);
        fmpz_mat_mul_classical_inline(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n\n");
            flint_printf("inputs:\n");
            fmpz_mat_print(A); flint_printf("\n\n");
            fmpz_mat_print(B); flint_printf("\n\n");
            flint_printf("fft:\n");
            fmpz_mat_print(C); flint_printf("\n\n");
            flint_printf("classical:\n");
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

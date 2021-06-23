/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mat.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("mul_fmpz_vec....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        fmpz * b, * c;
        slong j, m, n, blen;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        blen = n_randint(state, 50);

        fmpz_mat_init(C, m, 1);
        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, 1);
        c = _fmpz_vec_init(m);
        b = _fmpz_vec_init(blen);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(c, state, m, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(b, state, blen, n_randint(state, 200) + 1);

        fmpz_mat_mul_fmpz_vec(c, A, b, blen);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < n && j < blen; j++)
            fmpz_set(fmpz_mat_entry(B, j, 0), b + j);

        fmpz_mat_mul(C, A, B);

        for (j = 0; j < m; j++)
        {
            if (!fmpz_equal(fmpz_mat_entry(C, j, 0), c + j))
            {
                flint_printf("FAIL: wrong answer\n");
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        _fmpz_vec_clear(c, m);
        _fmpz_vec_clear(b, blen);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

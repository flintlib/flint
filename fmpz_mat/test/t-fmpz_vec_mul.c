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

    flint_printf("vec_fmpz_mul....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mat_t A, B, C;
        fmpz * a, * c;
        slong j, m, n, alen;

        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        fmpz_mat_init(C, 1, n);
        fmpz_mat_init(A, 1, m);
        fmpz_mat_init(B, m, n);
        c = _fmpz_vec_init(n);
        a = _fmpz_vec_init(alen);

        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(c, state, n, n_randint(state, 200) + 1);
        _fmpz_vec_randtest(a, state, alen, n_randint(state, 200) + 1);

        fmpz_mat_fmpz_vec_mul(c, a, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            fmpz_set(fmpz_mat_entry(A, 0, j), a + j);

        fmpz_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
        {
            if (!fmpz_equal(fmpz_mat_entry(C, 0, j), c + j))
            {
                flint_printf("FAIL: wrong answer\n");
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        _fmpz_vec_clear(c, n);
        _fmpz_vec_clear(a, alen);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

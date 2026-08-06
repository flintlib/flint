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
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fmpz_extras.h"

TEST_FUNCTION_START(fmpz_mat_charpoly_bound, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_t c, bound, height;
        fmpz_mat_t A;
        slong i, n, cbits, psteps, pbits;

        n = n_randint(state, 10);

        fmpz_init(c);
        fmpz_init(bound);
        fmpz_init(height);
        fmpz_mat_init(A, n, n);

        cbits = 1 + n_randint(state, 30);

        if (n_randint(state, 2) == 0)
        {
            fmpz_poly_t f;
            fmpz_poly_init(f);
            fmpz_mat_randtest(A, state, cbits);
            fmpz_mat_charpoly_berkowitz(f, A);
            fmpz_poly_height(height, f);
            fmpz_poly_clear(f);
        }
        else
        {
            /* Similiarity transform of random companion matrix */
            psteps = n_randint(state, 10);
            pbits = 1 + n_randint(state, 10);

            for (i = 0; i < n - 1; i++)
                fmpz_one(fmpz_mat_entry(A, i, i + 1));
            fmpz_one(height);
            for (i = 0; i < n; i++)
            {
                fmpz_randtest(fmpz_mat_entry(A, n - 1, i), state, cbits);
                fmpz_abs(c, fmpz_mat_entry(A, n - 1, i));
                fmpz_max(height, height, c);
            }

            for (i = 0; i < psteps; i++)
            {
               fmpz_randtest(c, state, pbits);
               fmpz_mat_similarity(A, n_randint(state, n), c);
            }
        }

        fmpz_mat_charpoly_bound(bound, A);

        if (fmpz_cmp(height, bound) > 0)
        {
            flint_printf("FAIL\n");
            flint_printf("A = \n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("height = %{fmpz}\n", height);
            flint_printf("bound = %{fmpz}\n", bound);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_clear(bound);
        fmpz_clear(height);
        fmpz_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}

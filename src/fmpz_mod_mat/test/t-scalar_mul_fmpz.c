/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod_mat.h"

TEST_FUNCTION_START(fmpz_mod_mat_scalar_mul_fmpz, state)
{
    slong m, n, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A, B, C, D;
        fmpz_t mod, c, c1;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(c);
        fmpz_init(c1);
        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        fmpz_randtest(c, state, 200);

        fmpz_mod_mat_init(A, m, n, mod);
        fmpz_mod_mat_init(B, m, n, mod);
        fmpz_mod_mat_init(C, m, n, mod);
        fmpz_mod_mat_init(D, m, n, mod);

        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);

        fmpz_mod_mat_scalar_mul_fmpz(C, A, c);
        fmpz_set(c1, c);
        fmpz_sub_ui(c1, c1, 1);
        fmpz_mod_mat_scalar_mul_fmpz(D, A, c1);

        /* c*A - (c-1)*A == A */
        fmpz_mod_mat_sub(D, C, D);

        if (!fmpz_mod_mat_equal(A, D))
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        /* Aliasing */
        fmpz_mod_mat_scalar_mul_fmpz(C, A, c);
        fmpz_mod_mat_scalar_mul_fmpz(A, A, c);

        if (!fmpz_mod_mat_equal(A, C))
        {
            flint_printf("FAIL\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
        fmpz_mod_mat_clear(D);
        fmpz_clear(mod);
        fmpz_clear(c1);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}

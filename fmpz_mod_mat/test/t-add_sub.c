/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

int
main(void)
{
    slong m, n, rep;
    FLINT_TEST_INIT(state);

    flint_printf("add/sub/neg....");
    fflush(stdout);

    

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A;
        fmpz_mod_mat_t B;
        fmpz_mod_mat_t C;
        fmpz_t mod;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        fmpz_mod_mat_init(A, m, n, mod);
        fmpz_mod_mat_init(B, m, n, mod);
        fmpz_mod_mat_init(C, m, n, mod);

        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);

        fmpz_mod_mat_neg(C, A);
        fmpz_mod_mat_add(A, A, B);
        fmpz_mod_mat_sub(A, A, B);
        fmpz_mod_mat_neg(A, A);

        if (!fmpz_mod_mat_equal(A, C))
        {
            flint_printf("FAIL: matrices not equal!\n");
            abort();
        }

        fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
        fmpz_mod_mat_clear(C);
	fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

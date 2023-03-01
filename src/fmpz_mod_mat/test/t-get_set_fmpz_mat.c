/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2021 William Hart

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

    flint_printf("get/set_fmpz_mat....");
    fflush(stdout);


    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mod_mat_t A;
        fmpz_mod_mat_t B;
        fmpz_t mod;
        fmpz_mat_t C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        fmpz_mod_mat_init(A, m, n, mod);
        fmpz_mod_mat_init(B, m, n, mod);

        fmpz_mat_init(C, m, n);

        fmpz_mod_mat_randtest(A, state);

        fmpz_mod_mat_get_fmpz_mat(C, A);
	fmpz_mod_mat_set_fmpz_mat(B, C);

        if (!fmpz_mod_mat_equal(A, B))
        {
            flint_printf("FAIL: matrices not equal!\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(C);
	
	fmpz_mod_mat_clear(A);
        fmpz_mod_mat_clear(B);
	fmpz_clear(mod);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

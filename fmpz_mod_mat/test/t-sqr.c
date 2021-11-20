/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_mat.h"
#include "ulong_extras.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("sqr....");
    fflush(stdout);

    

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mat_t A, B, C;
        slong n;
	fmpz_t mod;

        n = n_randint(state, 20);

        fmpz_init(mod);
        fmpz_randtest_not_zero(mod, state, 200);
        fmpz_abs(mod, mod);

        fmpz_mod_mat_init(A, n, n, mod);
        fmpz_mod_mat_init(B, n, n, mod);
        fmpz_mod_mat_init(C, n, n, mod);

        fmpz_mod_mat_randtest(A, state);
        fmpz_mod_mat_randtest(B, state);

        /* Make sure noise in the output is ok */
        fmpz_mod_mat_randtest(B, state);

        fmpz_mod_mat_sqr(B, A);
        fmpz_mod_mat_mul(C, A, A);

        if (!fmpz_mod_mat_equal(C, B))
        {
            flint_printf("FAIL: results not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_sqr(A, A);

        if (!fmpz_mod_mat_equal(A, B))
        {
            flint_printf("FAIL: aliasing failed\n");
            fflush(stdout);
            flint_abort();
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

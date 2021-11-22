/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"

int
main(void)
{
    int iter;

    FLINT_TEST_INIT(state);

    flint_printf("one/is_one....");
    fflush(stdout);    

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_poly_mat_t A;
        slong m, n;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_randtest(A, state, n_randint(state, 5));
        nmod_poly_mat_one(A);

        if (!nmod_poly_mat_is_one(A))
        {
            flint_printf("FAIL: expected matrix to be one\n");
            fflush(stdout);
            flint_abort();
        }

        if (m > 0 && n > 0)
        {
            m = n_randint(state, m);
            n = n_randint(state, n);

            if (m != n)
                nmod_poly_randtest_not_zero(nmod_poly_mat_entry(A, m, n),
                    state, 5);
            else
                do { nmod_poly_randtest_not_zero(nmod_poly_mat_entry(A, m, n),
                    state, 5); }
                while (nmod_poly_is_one(nmod_poly_mat_entry(A, m, n)));

            if (nmod_poly_mat_is_one(A))
            {
                flint_printf("FAIL: expected matrix not to be one\n");
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_poly_mat_clear(A);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

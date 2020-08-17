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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("neg....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_sub(C, A, B);

        nmod_mat_neg(B, B);
        nmod_mat_add(D, A, B);

        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL\n");
            abort();
        }

        nmod_mat_neg(C, B);
        nmod_mat_neg(B, B);

        if (!nmod_mat_equal(C, B))
        {
            flint_printf("FAIL\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

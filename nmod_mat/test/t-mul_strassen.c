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
    slong i;
    FLINT_TEST_INIT(state);
    

    flint_printf("mul_strassen....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t mod = n_randtest_not_zero(state);

        slong m, k, n;

        m = n_randint(state, 400);
        k = n_randint(state, 400);
        n = n_randint(state, 400);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul_classical(C, A, B);
        nmod_mat_mul_strassen(D, A, B);

        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
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

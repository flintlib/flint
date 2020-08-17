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
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, mod, rep;
    FLINT_TEST_INIT(state);
    

    flint_printf("scalar_mul....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t c;

        m = n_randint(state, 20);
        n = n_randint(state, 20);
        mod = n_randtest_not_zero(state);

        c = n_randint(state, mod);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_scalar_mul(C, A, c);
        nmod_mat_scalar_mul(D, A, nmod_sub(c, UWORD(1), A->mod));

        /* c*A - (c-1)*A == A */
        nmod_mat_sub(D, C, D);

        if (!nmod_mat_equal(A, D))
        {
            flint_printf("FAIL\n");
            abort();
        }

        /* Aliasing */
        nmod_mat_scalar_mul(C, A, c);
        nmod_mat_scalar_mul(A, A, c);

        if (!nmod_mat_equal(A, C))
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

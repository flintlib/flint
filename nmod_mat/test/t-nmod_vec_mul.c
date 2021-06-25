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
#include "nmod_mat.h"

int main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("nmod_vec_mul....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t p;
        nmod_mat_t A, B, C;
        mp_limb_t * a, * c;
        slong j, m, n, alen;

        p = n_randtest_not_zero(state);
        m = n_randint(state, 50);
        n = n_randint(state, 50);
        alen = n_randint(state, 50);

        nmod_mat_init(C, 1, n, p);
        nmod_mat_init(A, 1, m, p);
        nmod_mat_init(B, m, n, p);
        c = _nmod_vec_init(n);
        a = _nmod_vec_init(alen);

        nmod_mat_randtest(B, state);
        _nmod_vec_randtest(c, state, n, B->mod);
        _nmod_vec_randtest(a, state, alen, B->mod);

        nmod_mat_nmod_vec_mul(c, a, alen, B);

        /* supposed to match mul of the chopped or zero-extended a */
        for (j = 0; j < m && j < alen; j++)
            nmod_mat_entry(A, 0, j) = a[j];

        nmod_mat_mul(C, A, B);

        for (j = 0; j < n; j++)
        {
            if (nmod_mat_entry(C, 0, j) != c[j])
            {
                flint_printf("FAIL: wrong answer\n");
                flint_abort();
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        _nmod_vec_clear(c);
        _nmod_vec_clear(a);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}


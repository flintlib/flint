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

    flint_printf("mul_nmod_vec....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mp_limb_t p;
        nmod_mat_t A, B, C;
        mp_limb_t * b, * c;
        slong j, m, n, blen;

        p = n_randtest_not_zero(state);
        m = n_randint(state, 50);
        n = n_randint(state, 50);
        blen = n_randint(state, 50);

        nmod_mat_init(C, m, 1, p);
        nmod_mat_init(A, m, n, p);
        nmod_mat_init(B, n, 1, p);
        c = _nmod_vec_init(m);
        b = _nmod_vec_init(blen);

        nmod_mat_randtest(A, state);
        _nmod_vec_randtest(c, state, m, A->mod);
        _nmod_vec_randtest(b, state, blen, A->mod);

        nmod_mat_mul_nmod_vec(c, A, b, blen);

        /* supposed to match mul of the chopped or zero-extended b */
        for (j = 0; j < n && j < blen; j++)
            nmod_mat_entry(B, j, 0) = b[j];

        nmod_mat_mul(C, A, B);

        for (j = 0; j < m; j++)
        {
            if (nmod_mat_entry(C, j, 0) != c[j])
            {
                flint_printf("FAIL: wrong answer\n");
                flint_abort();
            }
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        _nmod_vec_clear(c);
        _nmod_vec_clear(b);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

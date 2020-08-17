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

void
nmod_mat_mul_check(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, k;

    mp_limb_t s0, s1, s2;
    mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s0 = s1 = s2 = UWORD(0);

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, A->rows[i][k], B->rows[k][j]);
                add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
            }

            NMOD_RED(s2, s2, C->mod);
            NMOD_RED3(s0, s2, s1, s0, C->mod);
            C->rows[i][j] = s0;
        }
    }
}

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);


    flint_printf("mul....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D;
        mp_limb_t mod;

        slong m, k, n;

        m = n_randint(state, 50);
        k = n_randint(state, 50);
        n = n_randint(state, 50);

        /* We want to generate matrices with many entries close to half
           or full limbs with high probability, to stress overflow handling */
        switch (n_randint(state, 3))
        {
            case 0:
                mod = n_randtest_not_zero(state);
                break;
            case 1:
                mod = UWORD_MAX/2 + 1 - n_randbits(state, 4);
                break;
            case 2:
            default:
                mod = UWORD_MAX - n_randbits(state, 4);
                break;
        }

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, k, mod);
        nmod_mat_init(C, m, k, mod);
        nmod_mat_init(D, m, k, mod);

        if (n_randint(state, 2))
            nmod_mat_randtest(A, state);
        else
            nmod_mat_randfull(A, state);

        if (n_randint(state, 2))
            nmod_mat_randtest(B, state);
        else
            nmod_mat_randfull(B, state);

        nmod_mat_randtest(C, state);  /* make sure noise in the output is ok */

        nmod_mat_mul(C, A, B);
        nmod_mat_mul_check(D, A, B);

        if (!nmod_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            abort();
        }

        if (n == k)
        {
            nmod_mat_mul(A, A, B);

            if (!nmod_mat_equal(A, C))
            {
                flint_printf("FAIL: aliasing failed\n");
                flint_abort();
            }
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

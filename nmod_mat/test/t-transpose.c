/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

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
    len_t m, n, mod, mod2, rep;
    flint_rand_t state;
    flint_randinit(state);

    printf("transpose....");
    fflush(stdout);

    /* Rectangular transpose, same modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, n, m, mod);
        nmod_mat_init(C, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_transpose(B, A);
        nmod_mat_transpose(C, B);

        if (!nmod_mat_equal(C, A))
        {
            printf("FAIL: C != A\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
    }

    /* Rectangular transpose, different modulus */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, AT, B, BT, AT2;

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        mod = n_randtest_not_zero(state);
        mod2 = n_randtest_not_zero(state);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(AT, n, m, mod);
        nmod_mat_init(B, m, n, mod2);
        nmod_mat_init(BT, n, m, mod2);
        nmod_mat_init(AT2, n, m, mod2);

        nmod_mat_randtest(A, state);
        nmod_mat_set(B, A);

        nmod_mat_transpose(AT, A);
        nmod_mat_transpose(BT, B);

        nmod_mat_set(AT2, AT);

        if (!nmod_mat_equal(BT, AT2))
        {
            printf("FAIL: AT != BT\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(AT);
        nmod_mat_clear(AT2);
        nmod_mat_clear(B);
        nmod_mat_clear(BT);
    }

    /* Self-transpose */
    for (rep = 0; rep < 100 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B;

        m = n_randint(state, 20);
        mod = n_randtest_not_zero(state);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(B, m, m, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_set(B, A);

        nmod_mat_transpose(B, B);
        nmod_mat_transpose(B, B);

        if (!nmod_mat_equal(B, A))
        {
            printf("FAIL: B != A\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}

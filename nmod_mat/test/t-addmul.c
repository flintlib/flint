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
    len_t i;
    flint_rand_t state;
    flint_randinit(state);

    printf("addmul....");
    fflush(stdout);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mat_t A, B, C, D, T, E;
        mp_limb_t mod = n_randtest_not_zero(state);

        len_t m, k, n;

        m = n_randint(state, 100);
        k = n_randint(state, 100);
        n = n_randint(state, 100);

        /* Force Strassen test */
        if (i < 5)
        {
            m += 300;
            k += 300;
            n += 300;
        }

        nmod_mat_init(A, m, k, mod);
        nmod_mat_init(B, k, n, mod);
        nmod_mat_init(C, m, n, mod);
        nmod_mat_init(D, m, n, mod);
        nmod_mat_init(T, m, n, mod);
        nmod_mat_init(E, m, n, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);
        nmod_mat_randtest(C, state);

        nmod_mat_addmul(D, C, A, B);

        nmod_mat_mul(T, A, B);
        nmod_mat_add(E, C, T);

        if (!nmod_mat_equal(D, E))
        {
            printf("FAIL: results not equal\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            nmod_mat_print_pretty(E);
            abort();
        }

        /* Check aliasing */
        nmod_mat_addmul(C, C, A, B);

        if (!nmod_mat_equal(C, E))
        {
            printf("FAIL: results not equal (aliasing)\n");
            nmod_mat_print_pretty(A);
            nmod_mat_print_pretty(B);
            nmod_mat_print_pretty(C);
            nmod_mat_print_pretty(D);
            nmod_mat_print_pretty(E);
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_mat_clear(E);
        nmod_mat_clear(T);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}

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
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"


int
main(void)
{
    nmod_mat_t A, x, b, Ax;
    slong i, m, r;
    int solved;
    mp_limb_t mod;
    FLINT_TEST_INIT(state);
    

    flint_printf("solve_vec....");
    fflush(stdout);

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        m = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(b, m, 1, mod);
        nmod_mat_init(x, m, 1, mod);
        nmod_mat_init(Ax, m, 1, mod);

        nmod_mat_randrank(A, state, m);
        nmod_mat_randtest(b, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = nmod_mat_solve_vec(x->entries, A, b->entries);
        nmod_mat_mul(Ax, A, x);

        if (!nmod_mat_equal(Ax, b) || !solved)
        {
            flint_printf("FAIL:\n");
            flint_printf("Ax != b!\n");
            flint_printf("A:\n");
            nmod_mat_print_pretty(A);
            flint_printf("b:\n");
            nmod_mat_print_pretty(b);
            flint_printf("x:\n");
            nmod_mat_print_pretty(x);
            flint_printf("Ax:\n");
            nmod_mat_print_pretty(Ax);
            flint_printf("\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(b);
        nmod_mat_clear(x);
        nmod_mat_clear(Ax);
    }

    /* Test singular systems */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        m = 1 + n_randint(state, 20);
        r = n_randint(state, m);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        nmod_mat_init(b, m, 1, mod);
        nmod_mat_init(x, m, 1, mod);
        nmod_mat_init(Ax, m, 1, mod);

        nmod_mat_randrank(A, state, r);
        nmod_mat_randtest(b, state);

        /* Dense */
        if (n_randint(state, 2))
            nmod_mat_randops(A, 1+n_randint(state, 1+m*m), state);

        solved = nmod_mat_solve_vec(x->entries, A, b->entries);

        if (solved)
        {
            flint_printf("FAIL:\n");
            flint_printf("singular system was 'solved'\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(b);
        nmod_mat_clear(x);
        nmod_mat_clear(Ax);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

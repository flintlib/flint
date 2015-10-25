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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, rep;
    ulong mod;
    FLINT_TEST_INIT(state);

    flint_printf("charpoly....");
    fflush(stdout);

    

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A, B, C, D;
        nmod_poly_t f, g;

        m = n_randint(state, 10);
        n = m;

        mod = n_randprime(state, 6, 0);

        nmod_mat_init(A, m, n, mod);
        nmod_mat_init(B, m, n, mod);
        nmod_mat_init(C, m, m, mod);
        nmod_mat_init(D, n, n, mod);
        nmod_poly_init(f, mod);
        nmod_poly_init(g, mod);

        nmod_mat_randtest(A, state);
        nmod_mat_randtest(B, state);

        nmod_mat_mul(C, A, B);
        nmod_mat_mul(D, B, A);

        nmod_mat_charpoly(f, C);
        nmod_mat_charpoly(g, D);

        if (!nmod_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), nmod_mat_print_pretty(A), flint_printf("\n");
            flint_printf("Matrix B:\n"), nmod_mat_print_pretty(B), flint_printf("\n");
            flint_printf("cp(AB) = "), nmod_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(BA) = "), nmod_poly_print_pretty(g, "X"), flint_printf("\n");
            abort();
        }

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

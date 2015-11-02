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
#include "fmpq_poly.h"
#include "fmpq_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, rep;
    FLINT_TEST_INIT(state);

    flint_printf("minpoly....");
    fflush(stdout);

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpq_mat_t A;
        fmpq_poly_t f, g, q, r;

        m = n_randint(state, 4);
        n = m;

        fmpq_mat_init(A, m, n);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(q);
        fmpq_poly_init(r);

        fmpq_mat_randtest(A, state, 10);
        
        fmpq_mat_charpoly(f, A);
        fmpq_mat_minpoly(g, A);

        fmpq_poly_divrem(q, r, f, g);

        if (!fmpq_poly_is_zero(r))
        {
            flint_printf("FAIL: minpoly(A) doesn't divide charpoly(A).\n");
            flint_printf("Matrix A:\n"), fmpq_mat_print(A), flint_printf("\n");
            flint_printf("cp(A) = "), fmpq_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("mp(A) = "), fmpq_poly_print_pretty(g, "X"), flint_printf("\n");
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(q);
        fmpq_poly_clear(r);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

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
#include "fmpz.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong m, n, rep;
    FLINT_TEST_INIT(state);

    flint_printf("charpoly....");
    fflush(stdout);

    

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpz_mat_t A;
        fmpz_poly_t f, g, q, r;

        m = n_randint(state, 4);
        n = m;

        fmpz_mat_init(A, m, n);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        
        fmpz_mat_randtest(A, state, 10);
        
        fmpz_mat_minpoly(f, A);
        fmpz_mat_charpoly(g, A);

        fmpz_poly_divrem(q, r, g, f);

        if (!fmpz_poly_is_zero(r))
        {
            flint_printf("FAIL: minpoly(A) doesn't divide charpoly(A).\n");
            flint_printf("Matrix A:\n"), fmpz_mat_print(A), flint_printf("\n");
            flint_printf("mp(A) = "), fmpz_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(A) = "), fmpz_poly_print_pretty(g, "X"), flint_printf("\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

    flint_printf("charpoly....");
    fflush(stdout);

    

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        fmpq_mat_t A, B, C, D;
        fmpq_poly_t f, g;

        m = n_randint(state, 4);
        n = m;

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, m, m);
        fmpq_mat_init(D, n, n);
        fmpq_poly_init(f);
        fmpq_poly_init(g);

        fmpq_mat_randtest(A, state, 10);
        fmpq_mat_randtest(B, state, 10);

        fmpq_mat_mul(C, A, B);
        fmpq_mat_mul(D, B, A);

        fmpq_mat_charpoly(f, C);
        fmpq_mat_charpoly(g, D);

        if (!fmpq_poly_equal(f, g))
        {
            flint_printf("FAIL: charpoly(AB) != charpoly(BA).\n");
            flint_printf("Matrix A:\n"), fmpq_mat_print(A), flint_printf("\n");
            flint_printf("Matrix B:\n"), fmpq_mat_print(B), flint_printf("\n");
            flint_printf("cp(AB) = "), fmpq_poly_print_pretty(f, "X"), flint_printf("\n");
            flint_printf("cp(BA) = "), fmpq_poly_print_pretty(g, "X"), flint_printf("\n");
            abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

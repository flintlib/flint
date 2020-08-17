/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "fmpz.h"

int
main(void)
{
    slong i;

    FLINT_TEST_INIT(state);

    flint_printf("sqr....");
    fflush(stdout);    

    /* Check evaluation homomorphism */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, C;
        nmod_mat_t a, c, d;
        mp_limb_t x, mod;
        slong m, deg;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(C, m, m, mod);

        nmod_mat_init(a, m, m, mod);
        nmod_mat_init(c, m, m, mod);
        nmod_mat_init(d, m, m, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(C, state, deg);  /* noise in output */

        nmod_poly_mat_sqr(C, A);

        x = n_randint(state, 0);

        nmod_poly_mat_evaluate_nmod(a, A, x);
        nmod_poly_mat_evaluate_nmod(d, C, x);
        nmod_mat_mul(c, a, a);

        if (!nmod_mat_equal(c, d))
        {
            flint_printf("FAIL:\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("C:\n");
            nmod_poly_mat_print(C, "x");
            flint_printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(C);

        nmod_mat_clear(a);
        nmod_mat_clear(c);
        nmod_mat_clear(d);
    }

    /* Check aliasing B and A */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B;
        slong m, deg;
        mp_limb_t mod;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 20);
        deg = 1 + n_randint(state, 10);

        nmod_poly_mat_init(A, m, m, mod);
        nmod_poly_mat_init(B, m, m, mod);

        nmod_poly_mat_randtest(A, state, deg);
        nmod_poly_mat_randtest(B, state, deg);  /* noise in output */

        nmod_poly_mat_sqr(B, A);
        nmod_poly_mat_sqr(A, A);

        if (!nmod_poly_mat_equal(B, A))
        {
            flint_printf("FAIL (aliasing):\n");
            flint_printf("A:\n");
            nmod_poly_mat_print(A, "x");
            flint_printf("B:\n");
            nmod_poly_mat_print(B, "x");
            flint_printf("\n");
            abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

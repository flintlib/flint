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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("trace....");
    fflush(stdout);

    

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_mat_t A, B, AB, BA;
        nmod_poly_t trab, trba;
        mp_limb_t mod;
        slong m, n;

        mod = n_randtest_prime(state, 0);
        m = n_randint(state, 10);
        n = n_randint(state, 10);

        nmod_poly_mat_init(A, m, n, mod);
        nmod_poly_mat_init(B, n, m, mod);
        nmod_poly_mat_init(AB, m, m, mod);
        nmod_poly_mat_init(BA, n, n, mod);

        nmod_poly_init(trab, mod);
        nmod_poly_init(trba, mod);

        nmod_poly_mat_randtest(A, state, 1 + n_randint(state, 10));
        nmod_poly_mat_randtest(B, state, 1 + n_randint(state, 10));

        nmod_poly_mat_mul(AB, A, B);
        nmod_poly_mat_mul(BA, B, A);

        nmod_poly_mat_trace(trab, AB);
        nmod_poly_mat_trace(trba, BA);

        if (!nmod_poly_equal(trab, trba))
        {
            flint_printf("FAIL:\n");
            nmod_poly_mat_print(A, "x"), flint_printf("\n");
            nmod_poly_mat_print(B, "x"), flint_printf("\n");
            nmod_poly_mat_print(AB, "x"), flint_printf("\n");
            nmod_poly_mat_print(BA, "x"), flint_printf("\n");
            flint_printf("tr(AB): "),  nmod_poly_print(trab),    flint_printf("\n");
            flint_printf("tr(BA): "),  nmod_poly_print(trba),    flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(AB);
        nmod_poly_mat_clear(BA);
        nmod_poly_clear(trab);
        nmod_poly_clear(trba);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

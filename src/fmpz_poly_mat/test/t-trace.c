/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly_mat.h"

TEST_FUNCTION_START(fmpz_poly_mat_trace, state)
{
    slong i;

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_mat_t A, B, AB, BA;
        fmpz_poly_t trab, trba;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpz_poly_mat_init(A, m, n);
        fmpz_poly_mat_init(B, n, m);
        fmpz_poly_mat_init(AB, m, m);
        fmpz_poly_mat_init(BA, n, n);

        fmpz_poly_init(trab);
        fmpz_poly_init(trba);

        fmpz_poly_mat_randtest(A, state, 1 + n_randint(state, 10),
            1 + n_randint(state, 100));
        fmpz_poly_mat_randtest(B, state, 1 + n_randint(state, 10),
            1 + n_randint(state, 100));

        fmpz_poly_mat_mul(AB, A, B);
        fmpz_poly_mat_mul(BA, B, A);

        fmpz_poly_mat_trace(trab, AB);
        fmpz_poly_mat_trace(trba, BA);

        if (!fmpz_poly_equal(trab, trba))
        {
            flint_printf("FAIL:\n");
            fmpz_poly_mat_print(A, "x"), flint_printf("\n");
            fmpz_poly_mat_print(B, "x"), flint_printf("\n");
            fmpz_poly_mat_print(AB, "x"), flint_printf("\n");
            fmpz_poly_mat_print(BA, "x"), flint_printf("\n");
            flint_printf("tr(AB): "),  fmpz_poly_print(trab),    flint_printf("\n");
            flint_printf("tr(BA): "),  fmpz_poly_print(trba),    flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_mat_clear(A);
        fmpz_poly_mat_clear(B);
        fmpz_poly_mat_clear(AB);
        fmpz_poly_mat_clear(BA);
        fmpz_poly_clear(trab);
        fmpz_poly_clear(trba);
    }

    TEST_FUNCTION_END(state);
}

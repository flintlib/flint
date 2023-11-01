/*
    Copyright (C) 2016, Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
TEST_FUNCTION_START(fmpz_poly_remove_content_2exp, state)
{
    int iter;

    for (iter = 0; iter < 1000; iter++)
    {
        fmpz_poly_t f, g1, g2, g3;

        fmpz_poly_init(f);
        fmpz_poly_init(g1);
        fmpz_poly_init(g2);
        fmpz_poly_init(g3);

        fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
        _fmpz_poly_remove_content_2exp(f->coeffs, f->length);

        fmpz_poly_scalar_mul_ui(g1, f, 2);
        _fmpz_poly_remove_content_2exp(g1->coeffs, g1->length);

        fmpz_poly_scalar_mul_ui(g2, f, 4);
        _fmpz_poly_remove_content_2exp(g2->coeffs, g2->length);

        fmpz_poly_scalar_mul_ui(g3, f, 256);
        _fmpz_poly_remove_content_2exp(g3->coeffs, g3->length);

        if ( ! (fmpz_poly_equal(f, g1) && fmpz_poly_equal(f, g2) && fmpz_poly_equal(f, g3)) )
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "); fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("g1 = "); fmpz_poly_print(g1); flint_printf("\n\n");
            flint_printf("g2 = "); fmpz_poly_print(g2); flint_printf("\n\n");
            flint_printf("g3 = "); fmpz_poly_print(g3); flint_printf("\n\n");
            flint_printf("ERROR \n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g1);
        fmpz_poly_clear(g2);
        fmpz_poly_clear(g3);
    }

    TEST_FUNCTION_END(state);
}

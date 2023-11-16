/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_get_coeff_ptr, state)
{
    int i, result;

    for (i = 0; i < 1000; i++)
    {
        arb_poly_t A;
        arb_t a;
        slong n = n_randint(state, 100);

        arb_poly_init(A);
        arb_poly_randtest(A, state, n_randint(state, 100), 100, 10);
        arb_init(a);

        arb_poly_get_coeff_arb(a, A, n);

        result = n < arb_poly_length(A) ?
                     arb_equal(a, arb_poly_get_coeff_ptr(A, n)) :
                     arb_poly_get_coeff_ptr(A, n) == NULL;
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A = "), arb_poly_printd(A, 10), flint_printf("\n\n");
            flint_printf("a = "), arb_print(a), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_abort();
        }

        arb_poly_clear(A);
        arb_clear(a);
    }

    TEST_FUNCTION_END(state);
}

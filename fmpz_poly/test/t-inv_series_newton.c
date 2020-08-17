/*
    Copyright (C) 2009 William Hart

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
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("inv_series_newton....");
    fflush(stdout);

    /* Check Q^{-1} * Q is congruent 1 mod t^n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c, one;
        slong n = n_randint(state, 80) + 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(one);

        fmpz_poly_randtest_not_zero(a, state, n_randint(state, 80) + 1, 100);
        fmpz_poly_set_coeff_si(a, 0, n_randint(state, 2) ? 1 : -1);

        fmpz_poly_set_ui(one, 1);

        fmpz_poly_inv_series_newton(b, a, n);
        fmpz_poly_mullow(c, a, b, n);

        result = (fmpz_poly_equal(c, one));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(one);
    }

    /* Verify bug fix for the case Q = -1 mod (x) */
    {
        fmpz_poly_t a, b, c, one;
        slong n = 1;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(one);

        fmpz_poly_set_si(a, -1);
        fmpz_poly_set_ui(one, 1);

        fmpz_poly_inv_series_newton(b, a, n);
        fmpz_poly_mullow(c, a, b, n);

        result = (fmpz_poly_equal(c, one));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_print(c), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(one);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

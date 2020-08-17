/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

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

    flint_printf("div_basecase....");
    fflush(stdout);

    

    /* Compare with full division, no aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r, q2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(q2);

        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, a->length + 1) + 1, 100);

        fmpz_poly_divrem_basecase(q, r, a, b);
        fmpz_poly_div_basecase(q2, a, b);

        result = (fmpz_poly_equal(q, q2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fmpz_poly_print(q2), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
        fmpz_poly_clear(q2);
    }

    /* Check q and a alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpz_poly_div_basecase(q, a, b);
        fmpz_poly_div_basecase(a, a, b);

        result = (fmpz_poly_equal(a, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }

    /* Check q and b alias */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_randtest(a, state, n_randint(state, 50), 100);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 50) + 1, 100);

        fmpz_poly_div_basecase(q, a, b);
        fmpz_poly_div_basecase(b, a, b);

        result = (fmpz_poly_equal(b, q));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(q), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
    }


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

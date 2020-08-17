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

    flint_printf("pseudo_rem....");
    fflush(stdout);

    

    /* Compare with divrem */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, q, r, r2;
        ulong d, d2;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(q);
        fmpz_poly_init(r);
        fmpz_poly_init(r2);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_divrem(q, r, &d, a, b);
        fmpz_poly_pseudo_rem(r2, &d2, a, b);

        result = (fmpz_poly_equal(r, r2) && d == d2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            fmpz_poly_print(r2), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(q);
        fmpz_poly_clear(r);
        fmpz_poly_clear(r2);
    }

    /* Check r and a alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;
        ulong d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_rem(r, &d, a, b);
        fmpz_poly_pseudo_rem(a, &d, a, b);

        result = (fmpz_poly_equal(a, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }

    /* Check r and b alias */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, r;
        ulong d;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(r);
        fmpz_poly_randtest(a, state, n_randint(state, 100), 50);
        fmpz_poly_randtest_not_zero(b, state, n_randint(state, 100) + 1, 50);

        fmpz_poly_pseudo_rem(r, &d, a, b);
        fmpz_poly_pseudo_rem(b, &d, a, b);

        result = (fmpz_poly_equal(b, r));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(r), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(r);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

/*
    Copyright (C) 2012 Fredrik Johansson

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
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("rfac_ui... ");
    fflush(stdout);

    

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, r;
        ulong a;

        fmpz_init(x);
        fmpz_init(r);

        fmpz_randtest(x, state, 1 + n_randint(state, 200));
        a = n_randint(state, 100);

        fmpz_rfac_ui(r, x, a);
        fmpz_rfac_ui(x, x, a);

        result = fmpz_equal(r, x);

        if (!result)
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("x: "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("a = %wu\n\n", a);
            flint_printf("r: "); fmpz_print(r); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(r);
    }

    /* Check rf(x,a) * rf(x+a,b) = rf(x,a+b) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t x, xa, r1, r2, r1r2, r3;
        ulong a, b;

        fmpz_init(x);
        fmpz_init(xa);
        fmpz_init(r1);
        fmpz_init(r2);
        fmpz_init(r1r2);
        fmpz_init(r3);

        fmpz_randtest(x, state, 1 + n_randint(state, 200));
        a = n_randint(state, 100);
        b = n_randint(state, 100);
        fmpz_add_ui(xa, x, a);

        fmpz_rfac_ui(r1, x, a);
        fmpz_rfac_ui(r2, xa, b);
        fmpz_rfac_ui(r3, x, a+b);

        fmpz_mul(r1r2, r1, r2);

        result = fmpz_equal(r1r2, r3);

        if (!result)
        {
            flint_printf("FAIL\n\n");
            flint_printf("x: "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("a = %wu, b = %wu\n\n", a, b);
            flint_printf("rf(x,a): "); fmpz_print(r1); flint_printf("\n\n");
            flint_printf("rf(x+a,b): "); fmpz_print(r2); flint_printf("\n\n");
            flint_printf("rf(x,a+b): "); fmpz_print(r3); flint_printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(xa);
        fmpz_clear(r1);
        fmpz_clear(r2);
        fmpz_clear(r1r2);
        fmpz_clear(r3);
    }

    
    flint_printf("PASS\n");
    FLINT_TEST_CLEANUP(state);
    return EXIT_SUCCESS;
}

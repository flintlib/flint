/*
    Copyright (C) 2011 Sebastian Pancratz

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

    flint_printf("divides....");
    fflush(stdout);

    /* Random values */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q, p;
        int divides;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);
	fmpz_init(q);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        divides = fmpz_divides(q, b, a);
        fmpz_mul(p, a, q);

        result = ((divides && fmpz_equal(p, b)) ||
	         (!divides && fmpz_is_zero(q) && !fmpz_equal(p, b)));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides = %d\n", divides);
	    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("p = "); fmpz_print(p); flint_printf("\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
	fmpz_clear(p);
        fmpz_clear(q);
    }

    /* Random b multiple of a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q, p;
        int divides;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(p);
	fmpz_init(q);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);
        fmpz_mul(b, a, b);

        divides = fmpz_divides(q, b, a);
        fmpz_mul(p, a, q);

        result = (divides && fmpz_equal(p, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides = %d\n", divides);
	    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("p = "); fmpz_print(p); flint_printf("\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
	fmpz_clear(p);
        fmpz_clear(q);
    }

    /* Check aliasing of q and a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q;
        int divides1, divides2;

        fmpz_init(a);
        fmpz_init(b);
	fmpz_init(q);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        divides1 = fmpz_divides(q, b, a);
        divides2 = fmpz_divides(a, b, a);

        result = (divides1 == divides2 && fmpz_equal(q, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
	    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(q);
    }

    /* Check aliasing of q and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, q;
        int divides1, divides2;

        fmpz_init(a);
        fmpz_init(b);
	fmpz_init(q);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        divides1 = fmpz_divides(q, b, a);
        divides2 = fmpz_divides(b, b, a);

        result = (divides1 == divides2 && fmpz_equal(q, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("divides1 = %d, divides2 = %d\n", divides1, divides2);
	    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("q = "); fmpz_print(q); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(q);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


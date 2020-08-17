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

    flint_printf("resultant....");
    fflush(stdout);

    

    /* Just one specific test */
    {
        fmpz_poly_t f, g;
        fmpz_t a, b;
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_poly_set_str(f, "11  -15 -2 -2 17 0 0 6 0 -5 1 -1");
        fmpz_poly_set_str(g, "9  2 1 1 1 1 1 0 -1 -2");
        fmpz_poly_resultant(a, f, g);
        fmpz_set_str(b, "-44081924855067", 10);

        result = (fmpz_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("res(f, h)  = "), fmpz_print(a), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    /* Check that R(fg, h) = R(f, h) R(g, h) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b, c, d;
        fmpz_poly_t f, g, h, p;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);
        fmpz_poly_init(p);
        fmpz_poly_randtest(f, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(g, state, n_randint(state, 50), 100);
        fmpz_poly_randtest(h, state, n_randint(state, 10), 100);

        fmpz_poly_resultant(a, f, h);
        fmpz_poly_resultant(b, g, h);
        fmpz_mul(c, a, b);
        fmpz_poly_mul(p, f, g);
        fmpz_poly_resultant(d, p, h);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_print(f), flint_printf("\n\n");
            flint_printf("g = "), fmpz_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpz_poly_print(h), flint_printf("\n\n");
            flint_printf("res(f, h)  = "), fmpz_print(a), flint_printf("\n\n");
            flint_printf("res(g, h)  = "), fmpz_print(b), flint_printf("\n\n");
            flint_printf("res(fg, h) = "), fmpz_print(d), flint_printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);
        fmpz_poly_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

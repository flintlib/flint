/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("compose_series....");
    fflush(stdout);

    

    /* Check aliasing of the first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series(f, g, h, n);
        fmpq_poly_compose_series(g, g, h, n);

        result = (fmpq_poly_equal(f, g));
        if (!result)
        {
            flint_printf("FAIL (aliasing 1):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(g), flint_printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check aliasing of the second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose_series(f, g, h, n);
        fmpq_poly_compose_series(h, g, h, n);

        result = (fmpq_poly_equal(f, h));
        if (!result)
        {
            flint_printf("FAIL (aliasing 2):\n");
            fmpq_poly_print(f), flint_printf("\n\n");
            fmpq_poly_print(h), flint_printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Compare with compose */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, h, s, t;
        slong n;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_randtest(g, state, n_randint(state, 40), 80);
        fmpq_poly_randtest(h, state, n_randint(state, 20), 50);
        fmpq_poly_set_coeff_ui(h, 0, 0);
        n = n_randint(state, 20);

        fmpq_poly_compose(s, g, h);
        fmpq_poly_truncate(s, n);
        fmpq_poly_compose_series(f, g, h, n);

        result = (fmpq_poly_equal(f, s));
        if (!result)
        {
            flint_printf("FAIL (comparison):\n");
            flint_printf("n = %wd\n", n);
            flint_printf("g = "), fmpq_poly_print(g), flint_printf("\n\n");
            flint_printf("h = "), fmpq_poly_print(h), flint_printf("\n\n");
            flint_printf("f = "), fmpq_poly_print(f), flint_printf("\n\n");
            flint_printf("s = "), fmpq_poly_print(s), flint_printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

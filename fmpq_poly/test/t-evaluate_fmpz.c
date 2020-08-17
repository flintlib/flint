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
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("evaluate_fmpz....");
    fflush(stdout);

    

    /* Check that (f+g)(a) = f(a) + g(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpq_poly_t f, g, h;
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(a);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, n_randint(state, 100));

        fmpq_poly_evaluate_fmpz(x, f, a);
        fmpq_poly_evaluate_fmpz(y, g, a);
        fmpq_add(x, x, y);
        fmpq_poly_add(h, f, g);
        fmpq_poly_evaluate_fmpz(y, h, a);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpq_poly_debug(f), flint_printf("\n");
            flint_printf("g = "), fmpq_poly_debug(g), flint_printf("\n");
            flint_printf("a = "), fmpz_print(a), flint_printf("\n");
            flint_printf("f(a) + g(a) = "), fmpq_print(x), flint_printf("\n\n");
            flint_printf("(f + g)(a)  = "), fmpq_print(y), flint_printf("\n\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(a);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
    }

    /* Check that (f*g)(a) = f(a) * g(a) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        fmpq_poly_t f, g;
        fmpq_t x, y;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(a);
        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(g, state, n_randint(state, 100), 200);
        fmpz_randtest(a, state, n_randint(state, 100));

        fmpq_poly_evaluate_fmpz(x, f, a);
        fmpq_poly_evaluate_fmpz(y, g, a);
        fmpq_mul(x, x, y);
        fmpq_poly_mul(f, f, g);
        fmpq_poly_evaluate_fmpz(y, f, a);

        result = (fmpq_equal(x, y));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(a), flint_printf("\n\n");
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(a);
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("evaluate_fmpq_poly....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        fmpq_poly_t f, g, h;
        qqbar_t x, fx, gx, hx, y;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);

        qqbar_init(x);
        qqbar_init(fx);
        qqbar_init(gx);
        qqbar_init(hx);
        qqbar_init(y);

        qqbar_randtest(x, state, 4, 10);
        qqbar_randtest(y, state, 4, 10);

        fmpq_poly_randtest(f, state, 8, 10);
        fmpq_poly_randtest(g, state, 5, 10);
        fmpq_poly_add(h, f, g);

        qqbar_evaluate_fmpq_poly(fx, f, x);
        qqbar_evaluate_fmpq_poly(gx, g, x);
        qqbar_evaluate_fmpq_poly(hx, h, x);
        qqbar_add(y, fx, gx);

        if (!qqbar_equal(y, hx))
        {
            flint_printf("FAIL!\n");
            flint_printf("f = "); fmpq_poly_print(f); flint_printf("\n\n");
            flint_printf("g = "); fmpq_poly_print(g); flint_printf("\n\n");
            flint_printf("h = "); fmpq_poly_print(h); flint_printf("\n\n");
            flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
            flint_printf("fx = "); qqbar_print(fx); flint_printf("\n\n");
            flint_printf("gx = "); qqbar_print(gx); flint_printf("\n\n");
            flint_printf("hx = "); qqbar_print(hx); flint_printf("\n\n");
            flint_printf("y = "); qqbar_print(y); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);

        qqbar_clear(x);
        qqbar_clear(fx);
        qqbar_clear(gx);
        qqbar_clear(hx);
        qqbar_clear(y);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


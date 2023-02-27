/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("shift_left/right....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 1000; i++)
    {
        arb_poly_t a, b;
        slong shift = n_randint(state, 100);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        arb_poly_shift_left(b, a, shift);
        arb_poly_shift_left(a, a, shift);

        result = (arb_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            arb_poly_printd(a, 10), flint_printf("\n\n");
            arb_poly_printd(b, 10), flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 1000; i++)
    {
        arb_poly_t a, b;
        slong shift = n_randint(state, 100);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        arb_poly_shift_right(b, a, shift);
        arb_poly_shift_right(a, a, shift);

        result = (arb_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            arb_poly_printd(a, 10), flint_printf("\n\n");
            arb_poly_printd(b, 10), flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 1000; i++)
    {
        arb_poly_t a, b, c;
        slong shift = n_randint(state, 100);

        arb_poly_init(a);
        arb_poly_init(b);
        arb_poly_init(c);
        arb_poly_randtest(a, state, n_randint(state, 100), 2 + n_randint(state, 200), 10);

        arb_poly_shift_left(b, a, shift);
        arb_poly_shift_right(c, b, shift);

        result = (arb_poly_equal(c, a));
        if (!result)
        {
            flint_printf("FAIL:\n");
            arb_poly_printd(a, 10), flint_printf("\n\n");
            arb_poly_printd(b, 10), flint_printf("\n\n");
            arb_poly_printd(c, 10), flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        arb_poly_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}


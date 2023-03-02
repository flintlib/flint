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
    ulong cflags = UWORD(0);

    FLINT_TEST_INIT(state);

    flint_printf("shift_left/right....");
    fflush(stdout);    

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong shift = (slong) n_randint(state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_shift_left(b, a, shift);
        fmpq_poly_shift_left(a, a, shift);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong shift;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 100) + 1, 200);

        shift = (slong) n_randint(state, a->length);

        fmpq_poly_shift_right(b, a, shift);
        fmpq_poly_shift_right(a, a, shift);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong shift = (slong) n_randint(state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_shift_left(b, a, shift);
        fmpq_poly_shift_right(c, b, shift);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(c, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

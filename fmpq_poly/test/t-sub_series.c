/*
    Copyright (C) 2009, 2014 William Hart

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

    flint_printf("sub_series....");
    fflush(stdout);    

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 200);

        fmpq_poly_sub_series(c, a, b, n);
        fmpq_poly_sub_series(a, a, b, n);
                cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
                result = (fmpq_poly_equal(a, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, c;
        slong n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(c);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 200);

        fmpq_poly_sub_series(c, a, b, n);
        fmpq_poly_sub_series(b, a, b, n);

        cflags |= fmpq_poly_is_canonical(b) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(c) ? 0 : 2;
        result = (fmpq_poly_equal(b, c) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(c), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(c);
    }

    /* Check truncate(a + b, n) = addseries(a, b, n) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, d, e;
        slong n;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(d);
        fmpq_poly_init(e);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);
        n = n_randint(state, 200);

        fmpq_poly_sub_series(d, a, b, n);
        fmpq_poly_sub(e, a, b);
        fmpq_poly_truncate(e, n);
        
        cflags |= fmpq_poly_is_canonical(d) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(e) ? 0 : 2;

        result = (fmpq_poly_equal(d, e) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fmpq_poly_debug(e), flint_printf("\n\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("cflags = %wu\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(d);
        fmpq_poly_clear(e);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

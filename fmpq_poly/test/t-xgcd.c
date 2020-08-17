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
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int cflags = 0, i, result;
    FLINT_TEST_INIT(state);
    
    flint_printf("xgcd....");
    fflush(stdout);

    

    /* Generic case, where a and b are coprime *******************************/

    /* Verify d == s a + t b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, d, e, f, s, t;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(d);
        fmpq_poly_init(e);
        fmpq_poly_init(f);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_randtest(a, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(b, state, n_randint(state, 60), 80);

        fmpq_poly_xgcd(d, s, t, a, b);
        fmpq_poly_mul(e, s, a);
        fmpq_poly_mul(f, t, b);
        fmpq_poly_add(e, e, f);

        cflags |= fmpq_poly_is_canonical(d) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(s) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(t) ? 0 : 4;
        result = (fmpq_poly_equal(d, e) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (correctness #1):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fmpq_poly_debug(s), flint_printf("\n\n");
            fmpq_poly_debug(t), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(d);
        fmpq_poly_clear(e);
        fmpq_poly_clear(f);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
    }

    /* Verify consistency with GCD */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, d, s, t;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(d);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_randtest(a, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(b, state, n_randint(state, 60), 80);

        fmpq_poly_xgcd(d, s, t, a, b);
        fmpq_poly_gcd(a, a, b);

        cflags |= fmpq_poly_is_canonical(d) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(s) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(t) ? 0 : 4;
        result = (fmpq_poly_equal(d, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (GCD #1):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fmpq_poly_debug(s), flint_printf("\n\n");
            fmpq_poly_debug(t), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(d);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
    }

    /* Generic case when a, b are most likely co-prime ***********************/

    /* Verify d == s a + t b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, d, s, t, z;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(d);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_init(z);
        fmpq_poly_randtest(a, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(b, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(z, state, n_randint(state, 20), 20);
        fmpq_poly_mul(a, a, z);
        fmpq_poly_mul(b, b, z);

        fmpq_poly_xgcd(d, s, t, a, b);
        fmpq_poly_mul(a, s, a);
        fmpq_poly_mul(b, t, b);
        fmpq_poly_add(a, a, b);

        cflags |= fmpq_poly_is_canonical(d) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(s) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(t) ? 0 : 4;
        result = (fmpq_poly_equal(d, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (correctness #2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fmpq_poly_debug(s), flint_printf("\n\n");
            fmpq_poly_debug(t), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(d);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
        fmpq_poly_clear(z);
    }

    /* Verify consistency with GCD */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, d, s, t, z;

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(d);
        fmpq_poly_init(s);
        fmpq_poly_init(t);
        fmpq_poly_init(z);
        fmpq_poly_randtest(a, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(b, state, n_randint(state, 60), 80);
        fmpq_poly_randtest(z, state, n_randint(state, 20), 20);
        fmpq_poly_mul(a, a, z);
        fmpq_poly_mul(b, b, z);

        fmpq_poly_xgcd(d, s, t, a, b);
        fmpq_poly_gcd(a, a, b);

        cflags |= fmpq_poly_is_canonical(d) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(s) ? 0 : 2;
        cflags |= fmpq_poly_is_canonical(t) ? 0 : 4;
        result = (fmpq_poly_equal(d, a) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (GCD #2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            fmpq_poly_debug(d), flint_printf("\n\n");
            fmpq_poly_debug(s), flint_printf("\n\n");
            fmpq_poly_debug(t), flint_printf("\n\n");
            flint_printf("cflags = %d\n\n", cflags);
            abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(d);
        fmpq_poly_clear(s);
        fmpq_poly_clear(t);
        fmpq_poly_clear(z);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


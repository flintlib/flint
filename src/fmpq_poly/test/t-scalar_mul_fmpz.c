/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "fmpz.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_scalar_mul_fmpz, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n;

        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), n_randint(state, 200));

        fmpq_poly_scalar_mul_fmpz(b, a, n);
        fmpq_poly_scalar_mul_fmpz(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that n (a + b) == na + nb */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpz_t n;

        fmpz_init(n);
        fmpz_randtest(n, state, 200);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), n_randint(state, 200));
        fmpq_poly_randtest(b, state, n_randint(state, 100), n_randint(state, 200));

        fmpq_poly_scalar_mul_fmpz(lhs, a, n);
        fmpq_poly_scalar_mul_fmpz(rhs, b, n);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_mul_fmpz(lhs, lhs, n);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Compare with fmpq_poly_scalar_mul_si */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpz_t n1;
        slong n;

        n = z_randtest(state);
        fmpz_init(n1);
        fmpz_set_si(n1, n);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), n_randint(state, 200));

        fmpq_poly_scalar_mul_fmpz(b, a, n1);
        fmpq_poly_scalar_mul_si(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n1);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}

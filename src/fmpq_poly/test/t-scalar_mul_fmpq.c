/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_scalar_mul_fmpq, state)
{
    int i, result;
    ulong cflags = UWORD(0);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        fmpq_t z;

        fmpq_init(z);
        fmpq_randtest(z, state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_mul_fmpq(b, a, z);
        fmpq_poly_scalar_mul_fmpq(a, a, z);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (aliasing):\n\n");
            flint_printf("a = "), fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("b = "), fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("z = "), fmpq_print(z), flint_printf("\n\n");
            gmp_printf("z = %Qd\n\n", z);
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check that (a * n1) * n2 == a * (n1 * n2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, lhs, rhs;
        fmpq_t z1, z2, z;

        fmpq_init(z1);
        fmpq_init(z2);
        fmpq_init(z);

        fmpq_randtest(z1, state, 100);
        fmpq_randtest(z2, state, 100);

        fmpq_mul(z, z1, z2);

        fmpq_poly_init(a);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_mul_fmpq(lhs, a, z1);
        fmpq_poly_scalar_mul_fmpq(lhs, lhs, z2);
        fmpq_poly_scalar_mul_fmpq(rhs, a, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            flint_printf("FAIL (a * n1 * n2):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            flint_printf("z1 = "), fmpq_print(z1), flint_printf("\n\n");
            flint_printf("z2 = "), fmpq_print(z2), flint_printf("\n\n");
            flint_printf("z = "), fmpq_print(z), flint_printf("\n\n");
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(z1);
        fmpq_clear(z2);
        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    /* Check that (a + b) * n == a*n + b*n */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b, lhs, rhs;
        fmpq_t z;

        fmpq_init(z);

        fmpq_randtest(z, state, 100);

        fmpq_poly_init(a);
        fmpq_poly_init(b);
        fmpq_poly_init(lhs);
        fmpq_poly_init(rhs);
        fmpq_poly_randtest(a, state, n_randint(state, 100), 200);
        fmpq_poly_randtest(b, state, n_randint(state, 100), 200);

        fmpq_poly_scalar_mul_fmpq(lhs, a, z);
        fmpq_poly_scalar_mul_fmpq(rhs, b, z);
        fmpq_poly_add(rhs, lhs, rhs);
        fmpq_poly_add(lhs, a, b);
        fmpq_poly_scalar_mul_fmpq(lhs, lhs, z);

        cflags |= fmpq_poly_is_canonical(lhs) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(rhs) ? 0 : 2;
        result = (fmpq_poly_equal(lhs, rhs) && !cflags);
        if (!result)
        {
            flint_printf("FAIL ((a + b) / n):\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("z = "), fmpq_print(z), flint_printf("\n\n");
            fmpq_poly_debug(lhs), flint_printf("\n\n");
            fmpq_poly_debug(rhs), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(z);
        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
        fmpq_poly_clear(lhs);
        fmpq_poly_clear(rhs);
    }

    TEST_FUNCTION_END(state);
}

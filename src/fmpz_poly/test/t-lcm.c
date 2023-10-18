/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_lcm, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);

        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_lcm(a, b, c);
        fmpz_poly_lcm(b, b, c);

        result = (fmpz_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and b):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t a, b, c;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);

        fmpz_poly_randtest(b, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(c, state, n_randint(state, 40), 80);

        fmpz_poly_lcm(a, b, c);
        fmpz_poly_lcm(c, b, c);

        result = (fmpz_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL (aliasing a and c):\n");
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
    }

    /* Check that GCD(f, g) LCM(f, g) == f g */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t f, g, gcd, lcm, lhs, rhs;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(gcd);
        fmpz_poly_init(lcm);
        fmpz_poly_init(lhs);
        fmpz_poly_init(rhs);

        fmpz_poly_randtest(f, state, n_randint(state, 40), 80);
        fmpz_poly_randtest(g, state, n_randint(state, 40), 80);

        fmpz_poly_gcd(gcd, f, g);
        fmpz_poly_lcm(lcm, f, g);
        fmpz_poly_mul(lhs, gcd, lcm);
        fmpz_poly_mul(rhs, f, g);
        if (!fmpz_poly_is_zero(rhs) && fmpz_sgn(fmpz_poly_lead(rhs)) < 0)
            fmpz_poly_neg(rhs, rhs);

        result = (fmpz_poly_equal(lhs, rhs));
        if (!result)
        {
            flint_printf("FAIL (GCD(f, g) * LCM(f, g) == f * g):\n");
            fmpz_poly_print(f), flint_printf("\n");
            fmpz_poly_print(g), flint_printf("\n");
            fmpz_poly_print(gcd), flint_printf("\n");
            fmpz_poly_print(lcm), flint_printf("\n");
            fmpz_poly_print(lhs), flint_printf("\n");
            fmpz_poly_print(rhs), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(gcd);
        fmpz_poly_clear(lcm);
        fmpz_poly_clear(lhs);
        fmpz_poly_clear(rhs);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz_poly.h"
#include "arb_poly.h"

TEST_FUNCTION_START(arb_poly_swinnerton_dyer_ui, state)
{
    slong iter;

    for (iter = 0; iter < 50 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t a, b;
        arb_t x, y;
        fmpz_poly_t c;
        slong i, n, prec;

        arb_poly_init(a);
        arb_poly_init(b);
        fmpz_poly_init(c);
        arb_init(x);
        arb_init(y);

        n = n_randint(state, 10);
        arb_poly_swinnerton_dyer_ui(a, n, 0);
        prec = 2 + n_randint(state, 10000);

        if (!arb_poly_get_unique_fmpz_poly(c, a))
        {
            flint_printf("FAIL (uniqueness)\n\n");
            flint_abort();
        }

        arb_poly_set_fmpz_poly(b, c, prec);

        arb_zero(x);
        for (i = 0; i < n; i++)
        {
            arb_sqrt_ui(y, n_nth_prime(i + 1), prec);
            if (n_randint(state, 2))
                arb_add(x, x, y, prec);
            else
                arb_sub(x, x, y, prec);
        }

        arb_poly_evaluate(y, b, x, prec);

        if (!arb_contains_zero(y))
        {
            flint_printf("FAIL (containment)\n\n");
            flint_abort();
        }

        arb_poly_swinnerton_dyer_ui(b, n, 2 + n_randint(state, 1000));

        if (!arb_poly_overlaps(a, b))
        {
            flint_printf("FAIL (overlap)\n\n");
            flint_abort();
        }

        arb_poly_clear(a);
        arb_poly_clear(b);
        fmpz_poly_clear(c);
        arb_clear(x);
        arb_clear(y);
    }

    TEST_FUNCTION_END(state);
}

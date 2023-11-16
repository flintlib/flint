/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_extras.h"
#include "fmpzi.h"

TEST_FUNCTION_START(fmpzi_gcd_euclidean_improved, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpzi_t x, y, d, xd, yd, g;

        fmpzi_init(x);
        fmpzi_init(y);
        fmpzi_init(d);
        fmpzi_init(xd);
        fmpzi_init(yd);
        fmpzi_init(g);

        fmpzi_randtest(x, state, 2 + n_randint(state, 100));
        fmpzi_randtest(y, state, 2 + n_randint(state, 100));
        fmpzi_randtest(xd, state, 2 + n_randint(state, 100));
        fmpzi_randtest(yd, state, 2 + n_randint(state, 100));

        fmpzi_gcd_euclidean_improved(d, x, y);

        if (fmpzi_is_zero(x) && fmpzi_is_zero(y))
        {
            fmpzi_zero(g);
        }
        else
        {
            if (!fmpzi_is_zero(d))
            {
                fmpzi_divrem(x, g, x, d);
                fmpzi_divrem(y, g, y, d);
            }

            fmpzi_mul_i_pow_si(x, x, n_randint(state, 4));
            fmpzi_mul_i_pow_si(y, y, n_randint(state, 4));

            fmpzi_randtest(d, state, 2 + n_randint(state, 100));
            fmpzi_mul(xd, x, d);
            fmpzi_mul(yd, y, d);

            fmpzi_gcd_euclidean_improved(g, xd, yd);

            fmpzi_canonicalise_unit(d, d);
        }

        if (!fmpzi_equal(g, d))
        {
            flint_printf("FAIL\n");
            flint_printf("x = "); fmpzi_print(x); printf("\n");
            flint_printf("y = "); fmpzi_print(y); printf("\n");
            flint_printf("d = "); fmpzi_print(d); printf("\n");
            flint_printf("xd = "); fmpzi_print(xd); printf("\n");
            flint_printf("yd = "); fmpzi_print(yd); printf("\n");
            flint_printf("g = "); fmpzi_print(g); printf("\n");
            flint_abort();
        }

        fmpzi_clear(x);
        fmpzi_clear(y);
        fmpzi_clear(d);
        fmpzi_clear(xd);
        fmpzi_clear(yd);
        fmpzi_clear(g);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_perfect_sqrt, state)
{
    slong i;

    for (i = 0; i < 5000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, s, s2;
        int square, square2;
        slong bits;

        fmpz_init(a); fmpz_init(s); fmpz_init(s2);

        bits = n_randint(state, 20) == 0 ? 400000 : 300;
        fmpz_randtest(a, state, bits);
        if (n_randint(state, 2))
            fmpz_mul(a, a, a);      /* square */
        if (n_randint(state, 4) == 0)
            fmpz_add_ui(a, a, n_randint(state, 3));

        square2 = fmpz_is_square(a);
        if (square2)
            fmpz_sqrt(s2, a);

        if (n_randint(state, 2))
        {
            fmpz_set(s, a);
            square = fmpz_perfect_sqrt(s, s);
        }
        else
            square = fmpz_perfect_sqrt(s, a);

        if (square != square2 || (square && !fmpz_equal(s, s2)))
            TEST_FUNCTION_FAIL("square = %d, square2 = %d\na = %{fmpz}\ns = %{fmpz}\n", square, square2, a, s);

        fmpz_clear(a); fmpz_clear(s); fmpz_clear(s2);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2025 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_factor.h"

TEST_FUNCTION_START(fmpz_factor_mul, state)
{
    slong i, b;
    fmpz_t x, y;
    fmpz_factor_t fx, fy, fz;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_factor_init(fx);
    fmpz_factor_init(fy);
    fmpz_factor_init(fz);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        b = n_randint(state, 1000000);
        if (n_randint(state, 2))
            b = -b;
        fmpz_set_si(x, b);

        b = n_randint(state, 1000000);
        if (n_randint(state, 2))
            b = -b;
        fmpz_set_si(y, b);

        fmpz_factor(fx, x);
        fmpz_factor(fy, y);

        fmpz_factor_mul(fz, fx, fy);

        fmpz_mul(x, x, y);
        fmpz_factor_expand(y, fz);

        if (fmpz_cmp(x, y) != 0)
        {
            flint_printf("FAIL (factor mul)\n");
            fmpz_print(x); flint_printf(" ");
            fmpz_print(y); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_factor_clear(fx);
    fmpz_factor_clear(fy);
    fmpz_factor_clear(fz);

    TEST_FUNCTION_END(state);
}

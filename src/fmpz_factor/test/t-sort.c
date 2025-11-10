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

TEST_FUNCTION_START(fmpz_factor_sort, state)
{
    slong i, j;
    fmpz_t x, y, z;
    fmpz_factor_t factor;

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_factor_init(factor);

        for(j = 0; j < 10 * flint_test_multiplier(); j++)
        {
            fmpz_randprime(x, state, 2 + j, 1);
            _fmpz_factor_append(factor, x, 1 + n_randint(state, 5));
        }
        fmpz_factor_expand(y, factor);
        fmpz_factor_sort(factor);
        fmpz_factor_expand(z, factor);

        // check order
        int ord = 1;
        for (j = 0; j < factor->num - 1; j++)
        {
            ord &= (fmpz_cmp(factor->p + j, factor->p + j + 1) < 0);
        }

        // check equality
        if (!fmpz_equal(y, z) || !ord)
        {
            flint_printf("FAIL (factor sort)\n");
            fmpz_factor_print(factor); flint_printf("\n");
            fmpz_print(y); flint_printf(" ");
            fmpz_print(z); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_factor_clear(factor);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);

    TEST_FUNCTION_END(state);
}
/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "long_extras.h"

TEST_FUNCTION_START(fmpz_get_si, state)
{
    int i, result;
    fmpz_t x;

    fmpz_init(x);

    fmpz_set_si(x, COEFF_MIN);
    if (COEFF_IS_MPZ(*x) || fmpz_get_si(x) != COEFF_MIN)
    {
        flint_printf("FAIL: COEFF_MIN");
        fflush(stdout);
        flint_abort();
    }

    fmpz_set_si(x, COEFF_MAX);
    if (COEFF_IS_MPZ(*x) || fmpz_get_si(x) != COEFF_MAX)
    {
        flint_printf("FAIL: COEFF_MIN");
        fflush(stdout);
        flint_abort();
    }

    fmpz_set_si(x, WORD_MIN);
    if (!COEFF_IS_MPZ(*x) || fmpz_get_si(x) != WORD_MIN)
    {
        flint_printf("FAIL: WORD_MIN");
        fflush(stdout);
        flint_abort();
    }

    fmpz_set_si(x, WORD_MIN);
    if (!COEFF_IS_MPZ(*x) || fmpz_get_si(x) != WORD_MIN)
    {
        flint_printf("FAIL: WORD_MAX");
        fflush(stdout);
        flint_abort();
    }

    fmpz_clear(x);

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        slong b, c;

        b = z_randtest(state);

        fmpz_init(a);

        fmpz_set_si(a, b);
        c = fmpz_get_si(a);

        result = (b == c) && _fmpz_is_canonical(a);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("b = %wd, c = %wd\n", b, c);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
    }

    TEST_FUNCTION_END(state);
}

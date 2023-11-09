/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_height, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        fmpz_t h;
        slong len, bits, bits2;

        fmpz_init(h);

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        bits = n_randint(state, 200);
        _fmpz_vec_randtest(a, state, len, bits);

        bits2 = _fmpz_vec_max_bits(a, len);
        _fmpz_vec_height(h, a, len);

        result = (fmpz_bits(h) == FLINT_ABS(bits2)) && (fmpz_sgn(h) >= 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wd, bits2 = %wd\n", bits, bits2);
            flint_printf("Computed height:\n");
            fmpz_print(h);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(h);
        _fmpz_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

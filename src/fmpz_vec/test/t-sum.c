/*
    Copyright (C) 2009, 2010 William Hart
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

TEST_FUNCTION_START(fmpz_vec_sum, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t x, y, z;

        slong len1 = n_randint(state, 100);
        slong len2 = n_randint(state, 100);

        a = _fmpz_vec_init(len1 + len2);
        b = a + len1;

        _fmpz_vec_randtest(a, state, len1 + len2, 200);

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);

        _fmpz_vec_sum(x, a, len1);
        _fmpz_vec_sum(y, b, len2);
        fmpz_add(x, x, y);
        _fmpz_vec_sum(z, a, len1 + len2);

        result = (fmpz_equal(x, z));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len1), flint_printf("\n\n");
            _fmpz_vec_print(b, len2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len1 + len2);

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
    }

    TEST_FUNCTION_END(state);
}

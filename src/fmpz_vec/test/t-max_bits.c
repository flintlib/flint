/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_max_bits, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len, bits, bits2, bits3;

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        bits = n_randint(state, 200);
        _fmpz_vec_randtest(a, state, len, bits);

        bits2 = _fmpz_vec_max_bits(a, len);
        bits3 = _fmpz_vec_max_bits_ref(a, len);

        result = (bits >= FLINT_ABS(bits2) && bits2 == bits3);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wd, bits2 = %wd bits3 = %wd\n", bits, bits2, bits3);
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

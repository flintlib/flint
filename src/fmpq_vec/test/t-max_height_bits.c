/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpq_vec_max_height_bits, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq *a;
        slong len;
        flint_bitcnt_t bits, bits2;

        len = n_randint(state, 100);

        a = _fmpq_vec_init(len);
        bits = n_randint(state, 200) + 1;
        _fmpq_vec_randtest(a, state, len, bits);

        bits2 = _fmpq_vec_max_height_bits(a, len);

        result = (bits >= FLINT_ABS(bits2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits = %wd, bits2 = %wd bits3 = %wd\n", bits, bits2);
            fflush(stdout);
            flint_abort();
        }

        _fmpq_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

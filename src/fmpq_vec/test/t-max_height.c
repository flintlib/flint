/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpq_vec_max_height, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq * a;
        fmpz_t h;
        slong len;
        flint_bitcnt_t bits, bits2;

        fmpz_init(h);

        len = n_randint(state, 100);

        a = _fmpq_vec_init(len);
        bits = n_randint(state, 200) + 1;
        _fmpq_vec_randtest(a, state, len, bits);

        bits2 = _fmpq_vec_max_height_bits(a, len);
        _fmpq_vec_max_height(h, a, len);

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
        _fmpq_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

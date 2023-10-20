/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_max_limbs, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len, bits;
        mp_size_t limbs, limbs2;

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        bits = n_randint(state, 200);
        limbs = (bits + FLINT_BITS - 1) / FLINT_BITS;
        _fmpz_vec_randtest(a, state, len, bits);

        limbs2 = _fmpz_vec_max_limbs(a, len);

        result = (limbs >= limbs2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("bits   = %wd\n", bits);
            flint_printf("limbs  = %wd\n", limbs);
            flint_printf("a      = {"), _fmpz_vec_print(a, len), flint_printf("}\n");
            flint_printf("limbs2 = %wd\n", limbs2);
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

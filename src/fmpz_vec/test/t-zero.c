/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_zero, state)
{
    int i, result;

    /* Check it's zero */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        slong len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_zero(a, len);

        result = (_fmpz_vec_is_zero(a, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_vec.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(d_vec_zero, state)
{
    int i;

    /* Check it's zero */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int result;
        double *a;
        slong len = n_randint(state, 100);

        a = _d_vec_init(len);
        _d_vec_randtest(a, state, len, 0, 0);

        _d_vec_zero(a, len);

        result = (_d_vec_is_zero(a, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _d_vec_clear(a);
    }

    TEST_FUNCTION_END(state);
}

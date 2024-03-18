/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "d_vec.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(d_vec_add, state)
{
    int i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        double *a, *b, *c;
        slong len = n_randint(state, 100);

        a = _d_vec_init(len);
        b = _d_vec_init(len);
        c = _d_vec_init(len);
        _d_vec_randtest(a, state, len, 0, 0);
        _d_vec_randtest(b, state, len, 0, 0);

        _d_vec_add(c, a, b, len);
        _d_vec_add(a, a, b, len);

        result = (_d_vec_equal(a, c, len));
        if (!result)
            TEST_FUNCTION_FAIL("");

        _d_vec_clear(a);
        _d_vec_clear(b);
        _d_vec_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        double *a, *b, *c;
        slong len = n_randint(state, 100);

        a = _d_vec_init(len);
        b = _d_vec_init(len);
        c = _d_vec_init(len);
        _d_vec_randtest(a, state, len, 0, 0);
        _d_vec_randtest(b, state, len, 0, 0);

        _d_vec_add(c, a, b, len);
        _d_vec_add(b, a, b, len);

        result = (_d_vec_equal(b, c, len));
        if (!result)
            TEST_FUNCTION_FAIL("");

        _d_vec_clear(a);
        _d_vec_clear(b);
        _d_vec_clear(c);
    }

    TEST_FUNCTION_END(state);
}

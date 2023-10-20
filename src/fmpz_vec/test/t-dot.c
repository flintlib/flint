/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_dot, state)
{
    int i, result;

    /* Check aliasing */

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t res1, res2;
        slong len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_set(b, a, len);

        fmpz_init(res1);
        fmpz_init(res2);

        _fmpz_vec_dot(res1, a, a, len);
        _fmpz_vec_dot(res2, a, b, len);

        result = fmpz_equal(res1, res2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(res1);
        fmpz_clear(res2);
    }

    /* Check commutative law */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t res1, res2;
        slong len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);

        fmpz_init(res1);
        fmpz_init(res2);

        _fmpz_vec_dot(res1, a, b, len);
        _fmpz_vec_dot(res2, b, a, len);

        result = fmpz_equal(res1, res2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(res1);
        fmpz_clear(res2);
    }

    TEST_FUNCTION_END(state);
}

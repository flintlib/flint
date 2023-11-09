/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpq_vec_dot, state)
{
    int i, result;

    /* Check commutative law */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq *a, *b;
        fmpq_t res1, res2;
        slong len = n_randint(state, 100);

        a = _fmpq_vec_init(len);
        b = _fmpq_vec_init(len);
        _fmpq_vec_randtest(a, state, len, 200);
        _fmpq_vec_randtest(b, state, len, 200);

        fmpq_init(res1);
        fmpq_init(res2);

        _fmpq_vec_dot(res1, a, b, len);
        _fmpq_vec_dot(res2, b, a, len);

        result = fmpq_equal(res1, res2);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpq_vec_clear(a, len);
        _fmpq_vec_clear(b, len);
		fmpq_clear(res1);
		fmpq_clear(res2);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpq_vec_get_set_fmpz_vec, state)
{
    int iter;

    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong i, n;
        fmpq * a, * b;
        fmpz * c;
        fmpz_t d;

        n = n_randint(state, 20);

        a = _fmpq_vec_init(n);
        b = _fmpq_vec_init(n);
        c = _fmpz_vec_init(n);
        fmpz_init(d);

        _fmpq_vec_randtest(a, state, n, 1 + n_randint(state, 200));
        _fmpq_vec_randtest(b, state, n, 1 + n_randint(state, 200));
        _fmpz_vec_randtest(c, state, n, 1 + n_randint(state, 200));

        _fmpq_vec_get_fmpz_vec_fmpz(c, d, a, n);
        _fmpq_vec_set_fmpz_vec(b, c, n);

        for (i = 0; i < n; i++)
        {
            fmpq_div_fmpz(b + i, b + i, d);
            if (!fmpq_equal(b + i, a + i))
            {
                flint_printf("FAIL: wrong answer\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpq_vec_clear(a, n);
        _fmpq_vec_clear(b, n);
        _fmpz_vec_clear(c, n);
        fmpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}

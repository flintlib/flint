/*
    Copyright (C) 2009, 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_mul_2exp, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        slong len = n_randint(state, 100);
        ulong exp = n_randint(state, 200);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mul_2exp(b, a, len, exp);
        _fmpz_vec_scalar_mul_2exp(a, a, len, exp);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("exp = %wu\n", exp);
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
    }

    /* Check aliasing of (a*2^e1)*2^e2 equals a*2^(e1+e2) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        slong len = n_randint(state, 100);
        ulong e1 = n_randint(state, 200);
        ulong e2 = n_randint(state, 200);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mul_2exp(b, a, len, e1);
        _fmpz_vec_scalar_mul_2exp(b, b, len, e2);
        _fmpz_vec_scalar_mul_2exp(a, a, len, e1 + e2);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("e1 = %wu, e2 = %wu\n", e1, e2);
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
    }

    TEST_FUNCTION_END(state);
}

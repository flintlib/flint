/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "long_extras.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_submul_si_2exp, state)
{
    int i, result;

    /* Compare with alternative method of computation */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c, *d;
        slong len = n_randint(state, 100), x;
        flint_bitcnt_t exp;

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);
        d = _fmpz_vec_init(len);

        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);
        _fmpz_vec_set(c, b, len);

        x = z_randtest(state);
        exp = n_randint(state, 200);

        _fmpz_vec_scalar_submul_si_2exp(b, a, len, x, exp);
        _fmpz_vec_scalar_mul_2exp(d, a, len, exp);
        _fmpz_vec_scalar_submul_si(c, d, len, x);

        result = (_fmpz_vec_equal(b, c, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("x = %wd, exp = %wu\n", x, exp);
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            _fmpz_vec_print(c, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);
        _fmpz_vec_clear(d, len);
    }

    TEST_FUNCTION_END(state);
}

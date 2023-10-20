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
#include "long_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_mul_si, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        slong len = n_randint(state, 100);
        slong n = z_randtest(state);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mul_si(b, a, len, n);
        _fmpz_vec_scalar_mul_si(a, a, len, n);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
    }

    /* Check agreement with _fmpz */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        slong len, n;
        fmpz_t x;

        len = n_randint(state, 100);
        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        n = z_randtest(state);
        fmpz_init(x);
        fmpz_set_si(x, n);

        _fmpz_vec_scalar_mul_fmpz(b, a, len, x);
        _fmpz_vec_scalar_mul_si(a, a, len, n);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            flint_printf("%li\n\n", n);
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}

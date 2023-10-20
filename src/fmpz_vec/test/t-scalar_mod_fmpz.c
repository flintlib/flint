/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_mod_fmpz, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz *a, *b;
        slong len = n_randint(state, 100);

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 1);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mod_fmpz(b, a, len, p);
        _fmpz_vec_scalar_mod_fmpz(a, a, len, p);

        result = (_fmpz_vec_equal(a, a, len));
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
        fmpz_clear(p);
    }

    /* Check the result is reduced */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz *a, *b;
        slong j, len = n_randint(state, 100);

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 1);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mod_fmpz(b, a, len, p);

        result = 1;
        for (j = 0; j < len; j++)
            result &= (fmpz_sgn(b + j) >= 0 && fmpz_cmp(b + j, p) < 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fmpz_print(p), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_content, state)
{
    int i, result;

    /* Check that content(a f) = abs(a) content(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d;
        fmpz *f;
        slong len = n_randint(state, 100);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);
        _fmpz_vec_randtest(f, state, len, 200);
        fmpz_randtest(a, state, 100);

        _fmpz_vec_content(c, f, len);
        _fmpz_vec_scalar_mul_fmpz(f, f, len, a);
        fmpz_abs(a, a);
        fmpz_mul(c, a, c);
        _fmpz_vec_content(d, f, len);

        result = (fmpz_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(c), flint_printf("\n\n");
            fmpz_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
    }

    TEST_FUNCTION_END(state);
}

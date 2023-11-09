/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpfr_vec.h"

TEST_FUNCTION_START(mpfr_vec_set_equal, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpfr_ptr a;
        slong len = n_randint(state, 100);

        a = _mpfr_vec_init(len, 200);
        _mpfr_vec_randtest(a, state, len);

        _mpfr_vec_set(a, a, len);

        result = (_mpfr_vec_equal(a, a, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _mpfr_vec_clear(a, len);
    }

    /* Compare copied vectors */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpfr_ptr a, b;
        slong len = n_randint(state, 100);

        a = _mpfr_vec_init(len, 200);
        b = _mpfr_vec_init(len, 200);
        _mpfr_vec_randtest(a, state, len);

        _mpfr_vec_set(b, a, len);

        result = (_mpfr_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _mpfr_vec_clear(a, len);
        _mpfr_vec_clear(b, len);
    }

    /* Compare unequal vectors */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpfr_ptr a, b;
        slong len = n_randint(state, 100) + 1;
        slong coeff;

        a = _mpfr_vec_init(len, 200);
        b = _mpfr_vec_init(len, 200);
        _mpfr_vec_randtest(a, state, len);

        _mpfr_vec_set(b, a, len);
        coeff = n_randint(state, len);
        mpfr_add_ui(b + coeff, b + coeff, 1, MPFR_RNDN);

        result = (!_mpfr_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            fflush(stdout);
            flint_abort();
        }

        _mpfr_vec_clear(a, len);
        _mpfr_vec_clear(b, len);
    }

    TEST_FUNCTION_END(state);
}

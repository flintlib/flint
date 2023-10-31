/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "mpf_vec.h"

TEST_FUNCTION_START(mpf_vec_set_fmpz_vec, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int result;
        fmpz *a;
        mpf *d1, *d2;
        slong bits, len;

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        d1 = _mpf_vec_init(len, 200);
        d2 = _mpf_vec_init(len, 200);

        bits = 1 + n_randint(state, 200);

        _fmpz_vec_randtest(a, state, len, bits);

        _mpf_vec_set_fmpz_vec(d1, a, len / 2);
        _mpf_vec_set_fmpz_vec(d1 + len / 2, a + len / 2, (len + 1) / 2);
        _mpf_vec_set_fmpz_vec(d2, a, len);

        result = (_mpf_vec_equal(d1, d2, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _mpf_vec_clear(d1, len);
        _mpf_vec_clear(d2, len);
    }

    TEST_FUNCTION_END(state);
}

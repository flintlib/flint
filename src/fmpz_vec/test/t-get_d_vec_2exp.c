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

#include <math.h>
#include "test_helpers.h"
#include "d_vec.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_get_d_vec_2exp, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int result;
        fmpz *a;
        double *d1, *d2;
        slong bits, j, len, l1, l2, l3;

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        d1 = _d_vec_init(len);
        d2 = _d_vec_init(len);

        bits = 1 + n_randint(state, 200);

        _fmpz_vec_randtest(a, state, len, bits);

        l1 = _fmpz_vec_get_d_vec_2exp(d1, a, len / 2);
        l2 = _fmpz_vec_get_d_vec_2exp(d1 + len / 2, a + len / 2,
                                      (len + 1) / 2);
        l3 = _fmpz_vec_get_d_vec_2exp(d2, a, len);

        if (l1 < l2)
            for (j = 0; j < len / 2; j++)
                d1[j] = ldexp(d1[j], l1 - l2);

        if (l2 < l1)
            for (j = len / 2; j < len; j++)
                d1[j] = ldexp(d1[j], l2 - l1);

        result = (l3 == FLINT_MAX(l1, l2) && _d_vec_equal(d1, d2, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(a, len);
        _d_vec_clear(d1);
        _d_vec_clear(d2);
    }

    TEST_FUNCTION_END(state);
}

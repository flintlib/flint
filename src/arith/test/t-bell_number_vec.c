/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "arith.h"

TEST_FUNCTION_START(arith_bell_number_vec, state)
{
    fmpz * b1;
    fmpz * b2;
    slong n;

    const slong maxn = 1000;


    b1 = _fmpz_vec_init(maxn);
    b2 = _fmpz_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        arith_bell_number_vec_recursive(b1, n);
        arith_bell_number_vec_multi_mod(b2, n);

        if (!_fmpz_vec_equal(b1, b2, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }
    }

    _fmpz_vec_clear(b1, maxn);
    _fmpz_vec_clear(b2, maxn);

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"

TEST_FUNCTION_START(arith_euler_number_zeta, state)
{
    fmpz * ress;
    fmpz_t res;
    slong n, N;


    N = 50;

    ress = _fmpz_vec_init(N);
    arith_euler_number_vec(ress, N);

    for (n = 0; n < N; n++)
    {
        fmpz_init(res);

        arith_euler_number(res, n);
        if (!fmpz_equal(res, ress + n))
        {
            flint_printf("FAIL: n = %wd\n", n);
            flint_printf("Value: "); fmpz_print(res); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(res);
    }

    _fmpz_vec_clear(ress, N);

    TEST_FUNCTION_END(state);
}

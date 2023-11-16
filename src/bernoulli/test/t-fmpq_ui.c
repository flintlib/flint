/*
    Copyright (C) 2011, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arf.h"
#include "arith.h"
#include "bernoulli.h"

TEST_FUNCTION_START(bernoulli_fmpq_ui, state)
{
    fmpz * num1;
    fmpz * den1;
    fmpz_t num2;
    fmpz_t den2;
    slong n, N;

    N = 4000 * FLINT_MIN(1.0, 0.1 * flint_test_multiplier());

    num1 = _fmpz_vec_init(N);
    den1 = _fmpz_vec_init(N);
    fmpz_init(num2);
    fmpz_init(den2);

    _arith_bernoulli_number_vec_multi_mod(num1, den1, N);

    for (n = 0; n < N; n++)
    {
        _bernoulli_fmpq_ui(num2, den2, n);

        if (!fmpz_equal(num1 + n, num2))
        {
            flint_printf("FAIL: n = %wd, numerator\n", n);
            flint_printf("vec:    "); fmpz_print(num1 + n); flint_printf("\n");
            flint_printf("single: "); fmpz_print(num2); flint_printf("\n");
            flint_abort();
        }

        if (!fmpz_equal(den1 + n, den2))
        {
            flint_printf("FAIL: n = %wd, denominator\n", n);
            flint_printf("vec:    "); fmpz_print(den1 + n); flint_printf("\n");
            flint_printf("single: "); fmpz_print(den2); flint_printf("\n");
            flint_abort();
        }
    }

    _fmpz_vec_clear(num1, N);
    _fmpz_vec_clear(den1, N);
    fmpz_clear(num2);
    fmpz_clear(den2);

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2021 Fredrik Johansson

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

TEST_FUNCTION_START(bernoulli_fmpq_ui_multi_mod, state)
{
    fmpz * num1;
    fmpz * den1;
    fmpz_t num2;
    fmpz_t den2;
    slong n, N;
    double alpha;

    N = 1500 * FLINT_MIN(1.0, 0.1 * flint_test_multiplier());

    num1 = _fmpz_vec_init(N);
    den1 = _fmpz_vec_init(N);
    fmpz_init(num2);
    fmpz_init(den2);

    _arith_bernoulli_number_vec_multi_mod(num1, den1, N);

    for (n = 0; n < N; n++)
    {
        if (n_randint(state, 2))
            alpha = -1.0;
        else
            alpha = n_randint(state, 11) / (double) 10;

        _bernoulli_fmpq_ui_multi_mod(num2, den2, n, alpha);

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

    {
        _bernoulli_fmpq_ui_multi_mod(num1, den1, 10000, 0.3);
        _bernoulli_fmpq_ui_multi_mod(num2, den2, 10000, 1.0);

        if (!fmpz_equal(num1, num2) || !fmpz_equal(den1, den2))
        {
            flint_printf("FAIL: n = 10000\n", n);
            flint_printf("num1    "); fmpz_print(num1); flint_printf("\n");
            flint_printf("num2    "); fmpz_print(num2); flint_printf("\n");
            flint_printf("den1    "); fmpz_print(den1); flint_printf("\n");
            flint_printf("den2    "); fmpz_print(den2); flint_printf("\n");
            flint_abort();
        }

        flint_set_num_threads(2);

        _bernoulli_fmpq_ui_multi_mod(num1, den1, 30000, -1.0);
        _bernoulli_fmpq_ui_multi_mod(num2, den2, 30000, 0.8);

        if (!fmpz_equal(num1, num2) || !fmpz_equal(den1, den2))
        {
            flint_printf("FAIL: n = 30000\n", n);
            flint_printf("num1    "); fmpz_print(num1); flint_printf("\n");
            flint_printf("num2    "); fmpz_print(num2); flint_printf("\n");
            flint_printf("den1    "); fmpz_print(den1); flint_printf("\n");
            flint_printf("den2    "); fmpz_print(den2); flint_printf("\n");
            flint_abort();
        }

        flint_set_num_threads(3);

        _bernoulli_fmpq_ui_multi_mod(num1, den1, 80000, -1.0);
        _bernoulli_fmpq_ui_multi_mod(num2, den2, 80000, 0.5);

        if (!fmpz_equal(num1, num2) || !fmpz_equal(den1, den2))
        {
            flint_printf("FAIL: n = 80000\n", n);
            flint_printf("num1    "); fmpz_print(num1); flint_printf("\n");
            flint_printf("num2    "); fmpz_print(num2); flint_printf("\n");
            flint_printf("den1    "); fmpz_print(den1); flint_printf("\n");
            flint_printf("den2    "); fmpz_print(den2); flint_printf("\n");
            flint_abort();
        }
    }

    _fmpz_vec_clear(num1, N);
    _fmpz_vec_clear(den1, N);
    fmpz_clear(num2);
    fmpz_clear(den2);

    TEST_FUNCTION_END(state);
}

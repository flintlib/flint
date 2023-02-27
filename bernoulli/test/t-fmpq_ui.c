/*
    Copyright (C) 2011, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint/fmpz_vec.h"
#include "flint/arith.h"
#include "bernoulli.h"

int main()
{
    fmpz * num1;
    fmpz * den1;
    fmpz_t num2;
    fmpz_t den2;
    slong n, N;

    flint_printf("fmpq_ui....");
    fflush(stdout);

    N = 4000 * FLINT_MIN(1.0, arb_test_multiplier());

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

    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}


/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"


int main()
{
    fmpz * num1;
    fmpz * num2;
    fmpz * num3;
    fmpz * den1;
    fmpz * den2;
    fmpz * den3;
    slong i, n, N;

    FLINT_TEST_INIT(state);

    flint_printf("bernoulli_number_vec....");
    fflush(stdout);

    N = 2000;

    num1 = _fmpz_vec_init(N);
    num2 = _fmpz_vec_init(N);
    num3 = _fmpz_vec_init(N);
    den1 = _fmpz_vec_init(N);
    den2 = _fmpz_vec_init(N);
    den3 = _fmpz_vec_init(N);

    for (n = 0; n < N; n += (n<100) ? 1 : n/3)
    {
        _arith_bernoulli_number_vec_recursive(num1, den1, n);
        _arith_bernoulli_number_vec_multi_mod(num2, den2, n);
        _arith_bernoulli_number_vec_zeta(num3, den3, n);

        for (i = 0; i < n; i++)
        {
            if (!fmpz_equal(num1 + i, num2 + i) ||
                !fmpz_equal(num1 + i, num3 + i))
            {
                flint_printf("FAIL: n = %wd, numerator of B_%wd\n", n, i);
                flint_printf("recursive: "); fmpz_print(num1 + i); flint_printf("\n");
                flint_printf("multi_mod: "); fmpz_print(num2 + i); flint_printf("\n");
                flint_printf("zeta:      "); fmpz_print(num3 + i); flint_printf("\n");
                abort();
            }

            if (!fmpz_equal(den1 + i, den2 + i) ||
                !fmpz_equal(den1 + i, den3 + i))
            {
                flint_printf("FAIL: n = %wd, denominator of B_%wd\n", n, i);
                flint_printf("recursive: "); fmpz_print(den1 + i); flint_printf("\n");
                flint_printf("multi_mod: "); fmpz_print(den2 + i); flint_printf("\n");
                flint_printf("zeta:      "); fmpz_print(den3 + i); flint_printf("\n");
                abort();
            }
        }
    }

    _fmpz_vec_clear(num1, N);
    _fmpz_vec_clear(num2, N);
    _fmpz_vec_clear(num3, N);
    _fmpz_vec_clear(den1, N);
    _fmpz_vec_clear(den2, N);
    _fmpz_vec_clear(den3, N);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

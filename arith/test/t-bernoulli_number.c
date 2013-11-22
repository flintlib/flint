/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

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
#include "fmpq.h"

int main()
{
    fmpz * num1;
    fmpz * den1;
    fmpz_t num2;
    fmpz_t den2;
    slong n, N;

    FLINT_TEST_INIT(state);

    flint_printf("bernoulli_number....");
    fflush(stdout);

    N = 4000;

    num1 = _fmpz_vec_init(N);
    den1 = _fmpz_vec_init(N);
    fmpz_init(num2);
    fmpz_init(den2);

    _arith_bernoulli_number_vec_multi_mod(num1, den1, N);

    for (n = 0; n < N; n++)
    {
        _arith_bernoulli_number(num2, den2, n);

        if (!fmpz_equal(num1 + n, num2))
        {
            flint_printf("FAIL: n = %wd, numerator\n", n);
            flint_printf("vec:    "); fmpz_print(num1 + n); flint_printf("\n");
            flint_printf("single: "); fmpz_print(num2); flint_printf("\n");
            abort();
        }

        if (!fmpz_equal(den1 + n, den2))
        {
            flint_printf("FAIL: n = %wd, denominator\n", n);
            flint_printf("vec:    "); fmpz_print(den1 + n); flint_printf("\n");
            flint_printf("single: "); fmpz_print(den2); flint_printf("\n");
            abort();
        }
    }

    /* Check non underscore versions */
    do
    {
        slong N = 100;
        fmpq * x;
        fmpq_t t;

        fmpq_init(t);
        x = flint_malloc(sizeof(fmpq) * N);

        for (n = 0; n < N; n++)
            fmpq_init(x + n);

        arith_bernoulli_number_vec(x, N);
        for (n = 0; n < N; n++)
        {
            arith_bernoulli_number(t, n);
            if (!fmpq_equal(x + n, t))
            {
                flint_printf("FAIL!: n = %wd\n", n);
                abort();
            }
        }

        for (n = 0; n < N; n++)
            fmpq_clear(x + n);
        flint_free(x);
        fmpq_clear(t);

    } while (0);

    _fmpz_vec_clear(num1, N);
    _fmpz_vec_clear(den1, N);
    fmpz_clear(num2);
    fmpz_clear(den2);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

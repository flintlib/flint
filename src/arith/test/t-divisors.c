/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "arith.h"

void arith_divisors_naive(fmpz_poly_t p, slong n)
{
    slong k;
    slong i = 0;

    n = FLINT_ABS(n);

    fmpz_poly_zero(p);
    for (k = 1; k <= n; k++)
    {
        if (n % k == 0)
        {
            fmpz_poly_set_coeff_si(p, i, k);
            i++;
        }
    }
}

TEST_FUNCTION_START(arith_divisors, state)
{
    fmpz_t t;
    fmpz_poly_t a, b;
    slong n;


    fmpz_init(t);
    fmpz_poly_init(a);
    fmpz_poly_init(b);

    for (n = -1000; n < 1000; n++)
    {
        fmpz_set_si(t, n);
        arith_divisors(a, t);
        arith_divisors_naive(b, n);
        if (!fmpz_poly_equal(a, b))
        {
            flint_printf("FAIL:\n");
            flint_printf("wrong value for n=%wd\n", n);
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_clear(t);
    fmpz_poly_clear(a);
    fmpz_poly_clear(b);

    TEST_FUNCTION_END(state);
}

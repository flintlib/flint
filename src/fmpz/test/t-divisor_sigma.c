/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpz.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "arith.h"

void fmpz_sigma_naive(fmpz_t x, ulong k, ulong n)
{
    slong i = 0;

    fmpz_t t;
    fmpz_poly_t p;
    fmpz_init(t);
    fmpz_poly_init(p);
    fmpz_set_ui(t, n);
    arith_divisors(p, t);

    fmpz_zero(x);
    for (i = 0; i < p->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(t, p, i);
        fmpz_pow_ui(t, t, k);
        fmpz_add(x, x, t);
    }

    fmpz_clear(t);
    fmpz_poly_clear(p);
}
TEST_FUNCTION_START(fmpz_divisor_sigma, state)
{
    fmpz_t m, a, b;
    slong n, k;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(m);

    for (n = 0; n < 5000; n++)
    {
        for (k = 0; k < 10; k++)
        {
            fmpz_set_ui(m, n);
            fmpz_divisor_sigma(a, k, m);
            fmpz_sigma_naive(b, k, n);
            if (!fmpz_equal(a, b) || !_fmpz_is_canonical(a))
            {
                flint_printf("FAIL:\n");
                flint_printf("wrong value for n=%wd, k=%wd\n", n, k);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(m);

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpz_poly_mulmid_KS, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        fmpz_poly_t a, b, c, d;
        slong nlo, nhi, bits;
        int aliasing, result;

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(c);
        fmpz_poly_init(d);

        bits = 1 + n_randint(state, 200);
        fmpz_poly_randtest(a, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(b, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(c, state, n_randint(state, 10), bits);
        fmpz_poly_randtest(d, state, n_randint(state, 10), bits);

        nlo = n_randint(state, 10);
        nhi = n_randint(state, 10);
        aliasing = n_randint(state, 5);

        if (aliasing == 3 || aliasing == 4)
            fmpz_poly_set(b, a);

        fmpz_poly_mul_classical(c, a, b);
        fmpz_poly_shift_right(c, c, nlo);
        fmpz_poly_truncate(c, FLINT_MAX(0, nhi - nlo));

        if (aliasing == 0)
        {
            fmpz_poly_mulmid_KS(d, a, b, nlo, nhi);
        }
        else if (aliasing == 1)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_KS(d, d, b, nlo, nhi);
        }
        else if (aliasing == 2)
        {
            fmpz_poly_set(d, b);
            fmpz_poly_mulmid_KS(d, a, d, nlo, nhi);
        }
        else if (aliasing == 3)
        {
            fmpz_poly_mulmid_KS(d, a, a, nlo, nhi);
        }
        else if (aliasing == 4)
        {
            fmpz_poly_set(d, a);
            fmpz_poly_mulmid_KS(d, d, d, nlo, nhi);
        }

        result = (fmpz_poly_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL: fmpz_poly_mulmid_KS\n");
            flint_printf("aliasing = %d, nlo = %wd, nhi = %wd\n", aliasing, nlo, nhi);
            fmpz_poly_print(a), flint_printf("\n\n");
            fmpz_poly_print(b), flint_printf("\n\n");
            fmpz_poly_print(c), flint_printf("\n\n");
            fmpz_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(c);
        fmpz_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}

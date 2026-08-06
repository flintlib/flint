/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_mulmid, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_poly_t a, b, c, d;
        ulong m = n_randtest_not_zero(state);
        slong nlo, nhi;
        int aliasing, result;
        slong N;

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(d, m);

        if (n_randint(state, 10) == 0)
            N = 500;
        else
            N = 50;

        nmod_poly_randtest(a, state, n_randint(state, N));
        nmod_poly_randtest(b, state, n_randint(state, N));
        nmod_poly_randtest(c, state, n_randint(state, N));
        nmod_poly_randtest(d, state, n_randint(state, N));

        nlo = n_randint(state, N);
        nhi = n_randint(state, N);
        aliasing = n_randint(state, 5);

        if (aliasing == 3 || aliasing == 4)
            nmod_poly_set(b, a);

        nmod_poly_mul(c, a, b);
        nmod_poly_shift_right(c, c, nlo);
        nmod_poly_truncate(c, FLINT_MAX(0, nhi - nlo));

        if (aliasing == 0)
        {
            nmod_poly_mulmid(d, a, b, nlo, nhi);
        }
        else if (aliasing == 1)
        {
            nmod_poly_set(d, a);
            nmod_poly_mulmid(d, d, b, nlo, nhi);
        }
        else if (aliasing == 2)
        {
            nmod_poly_set(d, b);
            nmod_poly_mulmid(d, a, d, nlo, nhi);
        }
        else if (aliasing == 3)
        {
            nmod_poly_mulmid(d, a, a, nlo, nhi);
        }
        else if (aliasing == 4)
        {
            nmod_poly_set(d, a);
            nmod_poly_mulmid(d, d, d, nlo, nhi);
        }

        result = (nmod_poly_equal(c, d));
        if (!result)
        {
            flint_printf("FAIL: nmod_poly_mulmid\n");
            flint_printf("aliasing = %d, nlo = %wd, nhi = %wd\n", aliasing, nlo, nhi);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            nmod_poly_print(d), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }

    TEST_FUNCTION_END(state);
}

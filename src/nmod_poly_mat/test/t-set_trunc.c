/*
    Copyright (C) 2023 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_poly_mat.h"

void test_with_dims(ulong m, ulong n, flint_rand_t state)
{
    int result;

    nmod_poly_mat_t a, b, c;
    ulong p;
    slong len;

    p = n_randtest_prime(state, 0);

    nmod_poly_mat_init(a, m, n, p);
    nmod_poly_mat_init(b, m, n, p);
    nmod_poly_mat_init(c, m, n, p);

    nmod_poly_mat_randtest(a, state, n_randint(state, 100));
    nmod_poly_mat_randtest(b, state, n_randint(state, 100));
    len = n_randint(state, 50);

    nmod_poly_mat_set_trunc(b, a, len);

    nmod_poly_t poly;
    nmod_poly_init(poly, p);
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            nmod_poly_set_trunc(poly, nmod_poly_mat_entry(a, i, j), len);
            result = nmod_poly_equal(poly, nmod_poly_mat_entry(b, i, j));
            if (!result)
            {
                flint_printf("FAIL:\n");
                nmod_poly_mat_print(a, "X"), flint_printf("\n\n");
                nmod_poly_mat_print(b, "X"), flint_printf("\n\n");
                flint_printf("truncation length %wd\n\n", len);
                fflush(stdout);
                flint_abort();
            }
        }
    }

    nmod_poly_mat_set(c, a);
    nmod_poly_mat_truncate(c, len);

    result = (nmod_poly_mat_equal(b, c));
    if (!result)
    {
        flint_printf("FAIL:\n");
        nmod_poly_mat_print(a, "X"), flint_printf("\n\n");
        nmod_poly_mat_print(b, "X"), flint_printf("\n\n");
        nmod_poly_mat_print(c, "X"), flint_printf("\n\n");
        flint_printf("truncation length %wd\n\n", len);
        fflush(stdout);
        flint_abort();
    }

    nmod_poly_mat_set_trunc(a, a, len);

    result = (nmod_poly_mat_equal(a, c));
    if (!result)
    {
        flint_printf("FAIL (aliasing):\n");
        nmod_poly_mat_print(a, "X"), flint_printf("\n\n");
        nmod_poly_mat_print(c, "X"), flint_printf("\n\n");
        flint_printf("truncation length %wd\n\n", len);
        fflush(stdout);
        flint_abort();
    }

    nmod_poly_mat_clear(a);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(c);
}

TEST_FUNCTION_START(nmod_poly_mat_set_trunc, state)
{
    for (int i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        test_with_dims(2,5,state);
        test_with_dims(3,3,state);
        test_with_dims(5,2,state);
    }

    TEST_FUNCTION_END(state);
}

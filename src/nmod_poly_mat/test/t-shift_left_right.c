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

/* Check a << shift >> shift == a */
void test_with_dimensions1(ulong rdim, ulong cdim, flint_rand_t state)
{
    int result;

    nmod_poly_mat_t a, b;
    mp_limb_t n = n_randtest_not_zero(state);
    slong shift = n_randint(state, 100);

    nmod_poly_mat_init(a, rdim, cdim, n);
    nmod_poly_mat_init(b, rdim, cdim, n);
    nmod_poly_mat_randtest(a, state, n_randint(state, 100));

    nmod_poly_mat_shift_left(b, a, shift);
    nmod_poly_mat_shift_right(b, b, shift);

    result = (nmod_poly_mat_equal(a, b));
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("shift = %wd, rdim = %ld, cdim = %ld, n = %wu\n",
                shift, rdim, cdim, n);
        nmod_poly_mat_print(a, "X"), flint_printf("\n\n");
        nmod_poly_mat_print(b, "X"), flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_poly_mat_clear(a);
    nmod_poly_mat_clear(b);
}

/* Check a << shift >> shift == a aliasing the other way */
void test_with_dimensions2(ulong rdim, ulong cdim, flint_rand_t state)
{
    int result;

    nmod_poly_mat_t a, b, c;
    mp_limb_t n = n_randtest_not_zero(state);
    slong shift = n_randint(state, 100);

    nmod_poly_mat_init(a, rdim, cdim, n);
    nmod_poly_mat_init(b, rdim, cdim, n);
    nmod_poly_mat_init(c, rdim, cdim, n);
    nmod_poly_mat_randtest(c, state, n_randint(state, 100));

    nmod_poly_mat_set(a, c);
    nmod_poly_mat_shift_left(c, c, shift);
    nmod_poly_mat_shift_right(b, c, shift);

    result = (nmod_poly_mat_equal(a, b));
    if (!result)
    {
        flint_printf("FAIL:\n");
        flint_printf("shift = %wd, rdim = %ld, cdim = %ld, n = %wu\n",
                shift, rdim, cdim, n);
        nmod_poly_mat_print(a, "X"), flint_printf("\n\n");
        nmod_poly_mat_print(b, "X"), flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_poly_mat_clear(a);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(c);
}

TEST_FUNCTION_START(nmod_poly_mat_shift_left_right, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        test_with_dimensions1(2, 5, state);
        test_with_dimensions1(3, 3, state);
        test_with_dimensions1(5, 2, state);
    }
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        test_with_dimensions2(2, 5, state);
        test_with_dimensions2(3, 3, state);
        test_with_dimensions2(5, 2, state);
    }

    TEST_FUNCTION_END(state);
}

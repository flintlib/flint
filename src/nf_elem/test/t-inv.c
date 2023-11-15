/*
    Copyright (C) 2014 William Hart
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_inv, state)
{
    int i, result;

    /* test a*^-1 = 1 */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, ainv, p1;

        nf_init_randtest(nf, state, 25, 100);

        nf_elem_init(a, nf);
        nf_elem_init(ainv, nf);
        nf_elem_init(p1, nf);

        do {
           nf_elem_randtest_not_zero(a, state, 100, nf);
        } while (!_nf_elem_invertible_check(a, nf));

        nf_elem_inv(ainv, a, nf);
        nf_elem_mul(p1, ainv, a, nf);

        result = (nf_elem_is_one(p1, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("ainv = "); nf_elem_print_pretty(ainv, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(ainv, nf);
        nf_elem_clear(p1, nf);

        nf_clear(nf);
    }

    /* test aliasing a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b;

        nf_init_randtest(nf, state, 25, 100);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        do {
           nf_elem_randtest_not_zero(b, state, 100, nf);
        } while (!_nf_elem_invertible_check(b, nf));

        nf_elem_inv(a, b, nf);
        nf_elem_inv(b, b, nf);

        result = (nf_elem_equal(a, b, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

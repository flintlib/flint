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

TEST_FUNCTION_START(nf_elem_div, state)
{
    int i, result;

    /* test a*^-1 = 1 */
    for (i = 0; i < 10*flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;

        nf_init_randtest(nf, state, 25, 100);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        do {
           nf_elem_randtest_not_zero(a, state, 100, nf);
        } while (!_nf_elem_invertible_check(a, nf));
        nf_elem_randtest(b, state, 100, nf);

        nf_elem_div(c, b, a, nf);
        nf_elem_mul(c, c, a, nf);

        result = (nf_elem_equal(b, c, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);
    }

    /* test aliasing a and b */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;

        nf_init_randtest(nf, state, 25, 100);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 100, nf);
        do {
           nf_elem_randtest_not_zero(c, state, 100, nf);
        } while (!_nf_elem_invertible_check(c, nf));

        nf_elem_div(a, b, c, nf);
        nf_elem_div(b, b, c, nf);

        result = (nf_elem_equal(a, b, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);
    }

    /* test aliasing a and c */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t pol;
        nf_t nf;
        nf_elem_t a, b, c;

        fmpq_poly_init(pol);
        do {
           fmpq_poly_randtest_not_zero(pol, state, 25, 100);
        } while (fmpq_poly_degree(pol) < 1);

        nf_init(nf, pol);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 100, nf);

        do {
           nf_elem_randtest_not_zero(c, state, 100, nf);
        } while (!_nf_elem_invertible_check(c, nf));

        nf_elem_div(a, b, c, nf);
        nf_elem_div(c, b, c, nf);

        result = (nf_elem_equal(a, c, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("d = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);

        fmpq_poly_clear(pol);
    }

    TEST_FUNCTION_END(state);
}

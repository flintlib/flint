/*
    Copyright (C) 2013 William Hart
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_mul, state)
{
    int i, result;

    /* test a*(b + c) = a*b + a*c */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c, s, p, p1, p2;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);
        nf_elem_init(s, nf);
        nf_elem_init(p, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);

        nf_elem_add(s, b, c, nf);
        nf_elem_mul(p1, a, b, nf);
        nf_elem_mul(p2, a, c, nf);

        nf_elem_mul(p, a, s, nf);
        nf_elem_add(s, p1, p2, nf);

        result = (nf_elem_equal(p, s, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("s = "); nf_elem_print_pretty(s, nf, "x"); printf("\n");
           printf("p = "); nf_elem_print_pretty(p, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);
        nf_elem_clear(s, nf);
        nf_elem_clear(p, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);

        nf_clear(nf);
    }

    /* test aliasing a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);

        nf_elem_mul(a, b, c, nf);
        nf_elem_mul(b, b, c, nf);

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
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);

        nf_elem_mul(a, b, c, nf);
        nf_elem_mul(c, b, c, nf);

        result = (nf_elem_equal(a, c, nf));
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

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2018 Tommy Hofmann
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_mul_gen, state)
{
    int i, result;

    /* test mul_gen(b) = a * b, where a is the generator */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, p1, p2;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);

        nf_elem_gen(a, nf);

        nf_elem_randtest(b, state, 200, nf);

        nf_elem_mul_gen(p1, b, nf);
        nf_elem_mul(p2, b, a, nf);

        result = (nf_elem_equal(p1, p2, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("K = "); nf_print(nf); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
           printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);

        nf_clear(nf);
    }

    /* test aliasing b and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t b, c;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        nf_elem_randtest(b, state, 200, nf);
        nf_elem_randtest(c, state, 200, nf);

        nf_elem_mul_gen(b, c, nf);
        nf_elem_mul_gen(c, c, nf);

        result = (nf_elem_equal(b, c, nf));
        if (!result)
        {
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

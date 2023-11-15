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

TEST_FUNCTION_START(nf_elem_set_equal, state)
{
    int i, result;

    /* set a = b, check a == b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);

        nf_elem_randtest(b, state, 200, nf);

        nf_elem_set(a, b, nf);

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

    /* test aliasing a and b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a;

        nf_init_randtest(nf, state, 40, 200);

        nf_elem_init(a, nf);

        nf_elem_randtest(a, state, 200, nf);

        nf_elem_set(a, a, nf);

        result = (nf_elem_equal(a, a, nf));
        if (!result)
        {
           printf("FAIL:\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           flint_abort();
        }

        nf_elem_clear(a, nf);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

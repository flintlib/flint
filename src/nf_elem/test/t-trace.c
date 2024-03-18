/*
    Copyright (C) 2014 William Hart
                  2020 Julian Rüth

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_trace, state)
{
    int i, result;

    /* test trace(a + b) = trace(a) + trace(b) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;
        fmpq_t atrace, btrace, ctrace, ctrace2;

        nf_init_randtest(nf, state, 25, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        fmpq_init(atrace);
        fmpq_init(btrace);
        fmpq_init(ctrace);
        fmpq_init(ctrace2);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);

        nf_elem_add(c, a, b, nf);
        nf_elem_trace(atrace, a, nf);
        nf_elem_trace(btrace, b, nf);
        nf_elem_trace(ctrace, c, nf);
        fmpq_add(ctrace2, atrace, btrace);

        result = (fmpq_equal(ctrace, ctrace2));
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf->pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("trace(a) = "); fmpq_print(atrace); printf("\n");
           printf("trace(b) = "); fmpq_print(btrace); printf("\n");
           printf("trace(a + b) = "); fmpq_print(ctrace); printf("\n");
           printf("trace(a) + trace(b) = "); fmpq_print(ctrace2); printf("\n");
           flint_abort();
        }

        fmpq_clear(atrace);
        fmpq_clear(btrace);
        fmpq_clear(ctrace);
        fmpq_clear(ctrace2);

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

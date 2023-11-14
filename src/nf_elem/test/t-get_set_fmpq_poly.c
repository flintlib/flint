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

TEST_FUNCTION_START(nf_elem_get_set_fmpq_poly, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f;
        fmpq_poly_t g;
        nf_t nf;
        nf_elem_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);

        nf_init_randtest(nf, state, 40, 200);

        fmpq_poly_randtest(f, state, fmpq_poly_degree(nf->pol) - 1, 200);

        nf_elem_init(a, nf);
        nf_elem_set_fmpq_poly(a, f, nf);
        nf_elem_get_fmpq_poly(g, a, nf);

        result = fmpq_poly_equal(f, g);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "a");
            printf("\n");
            flint_printf("f = "); fmpq_poly_print_pretty(f, "x");
            printf("\n");
            flint_printf("g = "); fmpq_poly_print_pretty(f, "x");
            flint_abort();
        }

        nf_elem_clear(a, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
    }

    /* try unreduced polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t f, g, r;
        nf_t nf;
        nf_elem_t a;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(r);

        nf_init_randtest(nf, state, 20, 200);

        fmpq_poly_randtest(f, state, n_randint(state, 30), 200);

        nf_elem_init(a, nf);
        nf_elem_set_fmpq_poly(a, f, nf);
        nf_elem_get_fmpq_poly(g, a, nf);

        fmpq_poly_rem(r, f, nf->pol);

        result = fmpq_poly_equal(r, g);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "a");
            printf("\n");
            flint_printf("f = "); fmpq_poly_print_pretty(f, "x");
            printf("\n");
            flint_printf("g = "); fmpq_poly_print_pretty(f, "x");
            printf("\n");
            flint_printf("r = "); fmpq_poly_print_pretty(r, "x");
            flint_abort();
        }

        nf_elem_clear(a, nf);

        nf_clear(nf);

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(r);
    }

    TEST_FUNCTION_END(state);
}

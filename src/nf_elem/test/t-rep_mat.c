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
#include "fmpq_mat.h"
#include "nf_elem.h"

TEST_FUNCTION_START(nf_elem_rep_mat, state)
{
    int i;

    /* test mul_gen(b) = a * b, where a is the generator */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, p1, p2, t;
        slong d;
        slong j, k;
        fmpq_mat_t R;

        nf_init_randtest(nf, state, 20, 100);

        d = fmpq_poly_degree(nf->pol);

        fmpq_mat_init(R, d, d);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(p1, nf);
        nf_elem_init(p2, nf);
        nf_elem_init(t, nf);

        nf_elem_randtest(b, state, 100, nf);

        nf_elem_rep_mat(R, b, nf);

        for (j = 0; j < d; j++)
        {
            nf_elem_gen(a, nf);
            nf_elem_pow(a, a, j, nf);
            nf_elem_mul(p1, b, a, nf);

            nf_elem_zero(p2, nf);

            for (k = 0; k < d; k++)
            {
                nf_elem_gen(t, nf);
                nf_elem_pow(t, t, k, nf);
                nf_elem_scalar_mul_fmpq(t, t, fmpq_mat_entry(R, j, k), nf);
                nf_elem_add(p2, p2, t, nf);
            }

            if (!nf_elem_equal(p1, p2, nf))
            {
                printf("FAIL:\n");
                printf("K = "); nf_print(nf); printf("\n");
                printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
                printf("p1 = "); nf_elem_print_pretty(p1, nf, "x"); printf("\n");
                printf("p2 = "); nf_elem_print_pretty(p2, nf, "x"); printf("\n");
                flint_abort();
            }
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(p1, nf);
        nf_elem_clear(p2, nf);
        nf_elem_clear(t, nf);
        fmpq_mat_clear(R);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

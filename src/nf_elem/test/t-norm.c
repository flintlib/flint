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

TEST_FUNCTION_START(nf_elem_norm, state)
{
    int i, result;

    /* test norm(a*b) = norm(a)*norm(b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b, c;
        fmpq_t anorm, bnorm, cnorm, cnorm2;

        nf_init_randtest(nf, state, 25, 200);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_init(c, nf);

        fmpq_init(anorm);
        fmpq_init(bnorm);
        fmpq_init(cnorm);
        fmpq_init(cnorm2);

        nf_elem_randtest(a, state, 200, nf);
        nf_elem_randtest(b, state, 200, nf);

        nf_elem_mul(c, a, b, nf);
        nf_elem_norm(anorm, a, nf);
        nf_elem_norm(bnorm, b, nf);
        nf_elem_norm(cnorm, c, nf);
        fmpq_mul(cnorm2, anorm, bnorm);

        result = (fmpq_equal(cnorm, cnorm2));
        if (!result)
        {
           printf("FAIL:\n");
           printf("nf->pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
           printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
           printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
           printf("c = "); nf_elem_print_pretty(c, nf, "x"); printf("\n");
           printf("norm(a) = "); fmpq_print(anorm); printf("\n");
           printf("norm(b) = "); fmpq_print(bnorm); printf("\n");
           printf("norm(a*b) = "); fmpq_print(cnorm); printf("\n");
           printf("norm(a)*norm(b) = "); fmpq_print(cnorm2); printf("\n");
           flint_abort();
        }

        fmpq_clear(anorm);
        fmpq_clear(bnorm);
        fmpq_clear(cnorm);
        fmpq_clear(cnorm2);

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        nf_elem_clear(c, nf);

        nf_clear(nf);
    }

    TEST_FUNCTION_END(state);
}

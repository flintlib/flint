/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"

TEST_FUNCTION_START(acb_poly_find_roots, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t A;
        acb_poly_t B;
        acb_poly_t C;
        acb_t t;
        acb_ptr roots;
        slong i, deg, isolated;
        slong prec = 10 + n_randint(state, 400);

        acb_init(t);
        acb_poly_init(A);
        acb_poly_init(B);
        acb_poly_init(C);

        do {
            acb_poly_randtest(A, state, 2 + n_randint(state, 15), prec, 5);
        } while (A->length == 0);
        deg = A->length - 1;

        roots = _acb_vec_init(deg);

        isolated = acb_poly_find_roots(roots, A, NULL, 0, prec);

        if (isolated == deg)
        {
            acb_poly_fit_length(B, 1);
            acb_set(B->coeffs, A->coeffs + deg);
            _acb_poly_set_length(B, 1);

            for (i = 0; i < deg; i++)
            {
                acb_poly_fit_length(C, 2);
                acb_one(C->coeffs + 1);
                acb_neg(C->coeffs + 0, roots + i);
                _acb_poly_set_length(C, 2);
                acb_poly_mul(B, B, C, prec);
            }

            if (!acb_poly_contains(B, A))
            {
                flint_printf("FAIL: product does not equal polynomial\n");
                acb_poly_printd(A, 15); flint_printf("\n\n");
                acb_poly_printd(B, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        for (i = 0; i < isolated; i++)
        {
            acb_poly_evaluate(t, A, roots + i, prec);
            if (!acb_contains_zero(t))
            {
                flint_printf("FAIL: poly(root) does not contain zero\n");
                acb_poly_printd(A, 15); flint_printf("\n\n");
                acb_printd(roots + i, 15); flint_printf("\n\n");
                acb_printd(t, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(roots, deg);

        acb_clear(t);
        acb_poly_clear(A);
        acb_poly_clear(B);
        acb_poly_clear(C);
    }

    TEST_FUNCTION_END(state);
}

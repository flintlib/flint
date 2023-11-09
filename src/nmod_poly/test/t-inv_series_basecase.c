/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_inv_series_basecase, state)
{
    int i, result;

    /* Check Q * Qinv = 1 mod x^n */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, qinv, prod;
        slong m;

        mp_limb_t n;
        do n = n_randtest_not_zero(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(prod, n);
        nmod_poly_init(qinv, n);
        nmod_poly_init(q, n);

        do nmod_poly_randtest(q, state, n_randint(state, 2000));
        while (q->length == 0 || q->coeffs[0] == 0);

        m = n_randint(state, 2000) + 1;

        nmod_poly_inv_series_basecase(qinv, q, m);

        nmod_poly_mul(prod, q, qinv);
        nmod_poly_truncate(prod, m);

        result = (prod->length == 1 && prod->coeffs[0] == 1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(qinv), flint_printf("\n\n");
            nmod_poly_print(prod), flint_printf("\n\n");
            flint_printf("n = %wd\n", n);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(qinv);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of q and qinv */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t q, qinv;
        slong m;

        mp_limb_t n;
        do n = n_randtest(state);
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(qinv, n);
        do nmod_poly_randtest(q, state, n_randint(state, 1000));
        while (q->length == 0 || q->coeffs[0] == 0);

        m = n_randint(state, 2000) + 1;

        nmod_poly_inv_series_basecase(qinv, q, m);
        nmod_poly_inv_series_basecase(q, q, m);

        result = (nmod_poly_equal(q, qinv));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(q), flint_printf("\n\n");
            nmod_poly_print(qinv), flint_printf("\n\n");
            nmod_poly_print(q), flint_printf("\n\n");
            flint_printf("n = %wd, m = %wd\n", n, m);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(qinv);
    }

    TEST_FUNCTION_END(state);
}

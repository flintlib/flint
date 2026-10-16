/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_num_real_roots_vca, state)
{
    slong iter;

    /* call with random nonzero polynomials */
    for (iter = 0; iter < 100; iter++)
    {
        slong k;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        fmpz_poly_randtest_not_zero(p, state, 20, 10 + n_randint(state, 100));
        k = fmpz_poly_num_real_roots_vca(p);
        if (k < 0 || k > fmpz_poly_degree(p))
        {
            flint_printf("ERROR:\n");
            flint_printf("got k in wrong range k = %wd\n", k);
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    /* we check on products of the form               */
    /*   Prod (X - r_i) x  R                          */
    /* where r_i are rationals and R has no real root */
    for (iter = 0; iter < 200; iter++)
    {
        slong k, n;
        fmpz_poly_t p,q;
        fmpq * vec;

        n = 1 + (slong)n_randint(state, 10);

        vec = _fmpq_vec_init(n);
        _fmpq_vec_randtest_uniq_sorted(vec, state, n, 80);

        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_product_roots_fmpq_vec(p, vec, n);
        fmpz_poly_randtest_no_real_root(q, state, 1 + (slong)n_randint(state, 5), 80);
        /* note: should we check that q is squarefree? */
        fmpz_poly_mul(p, p, q);

        k = fmpz_poly_num_real_roots_vca(p);
        if (k != n)
        {
            flint_printf("ERROR:\n");
            flint_printf("found k = %wd instead of n = %wd\n", k, n);
            flint_printf("vec = "); _fmpq_vec_print(vec, n); flint_printf("\n");
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_abort();
        }

        _fmpq_vec_clear(vec, n);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}

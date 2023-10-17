/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpz_poly_num_real_roots_sturm, state)
{
    int iter;

    /* call with random nonzero polynomials */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong k;
        fmpz_poly_t p;

        fmpz_poly_init(p);
        fmpz_poly_randtest_not_zero(p, state, 20, 10 + n_randint(state, 100));
        k = fmpz_poly_num_real_roots_sturm(p);
        if (k < 0 || k > fmpz_poly_degree(p))
        {
            printf("ERROR:\n");
            flint_printf("got k in wrong range k = %wd\n", k);
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(p);
    }

    /* we check on products of the form               */
    /*   Prod (X - r_i) x  R                          */
    /* where r_i are rationals and R has no real root */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong k, k_neg, k_pos, n;
        fmpz_poly_t p,q;
        fmpq * vec;

        n = 1 + n_randint(state, 10);

        vec = _fmpq_vec_init(n);
        _fmpq_vec_randtest_uniq_sorted(vec, state, n, 80);

        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_product_roots_fmpq_vec(p, vec, n);
        fmpz_poly_randtest_no_real_root(q, state, 1 + n_randint(state, 5), 80);
        /* note: here there is no need to check that q is squarefree (Sturm test remains valid) */
        fmpz_poly_mul(p, p, q);

        k = fmpz_poly_num_real_roots_sturm(p);
        if (k != n)
        {
            printf("ERROR:\n");
            flint_printf("found k = %wd instead of n = %wd\n", k, n);
            printf("vec = "); _fmpq_vec_print(vec, n); printf("\n");
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_is_zero(p->coeffs))
        {
            _fmpz_poly_num_real_roots_sturm(&k_neg, &k_pos, p->coeffs, p->length);
            for (k = 0; (k < n) && (fmpq_sgn(vec + k) < 0); k++);
            if ((k_neg + k_pos != n) || (k != k_neg))
            {
                printf("ERROR:\n");
                flint_printf("found k_neg = %wd and k_pos = %wd instead of %wd and %wd\n",
                        k_neg, k_pos, k, n - k);
                printf("vec = "); _fmpq_vec_print(vec, n); printf("\n");
                printf("p = "); fmpz_poly_print(p); printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpq_vec_clear(vec, n);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}

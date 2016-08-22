/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpq_vec.h"

int main()
{
    int iter;

    FLINT_TEST_INIT(state);

    printf("num_real_roots_sturm....");

    /* we check on products of the form               */
    /*   Prod (X - r_i) x  R                          */
    /* where r_i are rationals and R has no real root */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong k, n;
        fmpz_poly_t p,q;
        fmpq * vec;

        n = 1 + n_randint(state, 20);

        vec = _fmpq_vec_init(n);
        _fmpq_vec_randtest_uniq_sorted(vec, state, n, 30);

        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_set_rational_roots(p, vec, n);
        fmpz_poly_randtest_no_real_root(q, state, 1 + n_randint(state, 10), 50);
        /* note: here there is no need to check that q is squarefree (Sturm remains valid) */
        fmpz_poly_mul(p, p, q);

        k = fmpz_poly_num_real_roots_sturm(p);
        if (k != n)
        {
            printf("ERROR:\n");
            flint_printf("found k = %wd instead of n = %wd\n", k, n);
            printf("vec = "); _fmpq_vec_print(vec, n); printf("\n");
            printf("p = "); fmpz_poly_print(p); printf("\n");
            abort();
        }

        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

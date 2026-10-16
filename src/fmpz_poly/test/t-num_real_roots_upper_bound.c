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

TEST_FUNCTION_START(fmpz_poly_num_real_roots_upper_bound, state)
{
    slong iter;

    /* test polynomials with random rational roots */
    for( iter = 0; iter <= 1000; iter++ )
    {
        slong n_real_roots, n_complex_roots;
        fmpq * real_roots;

        slong bound;
        fmpz_poly_t p,q;

        n_real_roots = (slong)n_randint(state, 30);
        n_complex_roots = 1 + (slong)n_randint(state, 20);

        real_roots = _fmpq_vec_init(n_real_roots);

        _fmpq_vec_randtest(real_roots, state, n_real_roots, 100);

        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpz_poly_randtest_no_real_root(p, state, n_complex_roots, 40);
        fmpz_poly_product_roots_fmpq_vec(q, real_roots, n_real_roots);
        fmpz_poly_mul(p, p, q);

        bound = fmpz_poly_num_real_roots_upper_bound(p);

        if (n_real_roots > bound)
        {
            flint_printf("FAIL:\n");
            flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
            flint_printf("n_real_roots = %ld\n", n_real_roots);
            flint_printf("n_complex_roots  = %ld\n", n_complex_roots);
            flint_printf("got bound = %wd\n", bound);
            flint_abort();
        }

        _fmpq_vec_clear(real_roots, n_real_roots);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }

    TEST_FUNCTION_END(state);
}

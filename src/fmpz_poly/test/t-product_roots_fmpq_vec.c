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

TEST_FUNCTION_START(fmpz_poly_product_roots_fmpq_vec, state)
{
    slong iter;

    for (iter = 0; iter < 100; iter++)
    {
        fmpq * vec;
        fmpz_poly_t p;
        fmpq_t res;
        slong i, n;

        n = n_randint(state, 100);
        vec = _fmpq_vec_init(n);
        _fmpq_vec_randtest(vec, state, n, 100);

        fmpz_poly_init(p);
        fmpz_poly_product_roots_fmpq_vec(p, vec, n);

        if (fmpz_poly_degree(p) != n)
        {
            printf("ERROR:\n");
            flint_printf("expected degree %wd and got %wd",
                    n, fmpz_poly_degree(p));
            printf("vec = "); _fmpq_vec_print(vec, n); printf("\n");
            printf("p = "); fmpz_poly_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_init(res);

        for (i = 0; i < n; i++)
        {
            fmpz_poly_evaluate_fmpq(res, p, vec + i);
            if (!fmpq_is_zero(res))
            {
                printf("ERROR:\n");
                flint_printf("evaluation at the %wd-th root ", i);
                fmpq_print(vec + i);
                printf(" is not zero");
                printf("vec = "); _fmpq_vec_print(vec, n); printf("\n");
                printf("p = "); fmpz_poly_print(p); printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_clear(res);
        _fmpq_vec_clear(vec, n);
        fmpz_poly_clear(p);
    }

    TEST_FUNCTION_END(state);
}

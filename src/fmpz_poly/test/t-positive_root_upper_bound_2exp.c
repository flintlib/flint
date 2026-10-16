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

TEST_FUNCTION_START(fmpz_poly_positive_root_upper_bound_2exp, state)
{
    slong iter;

    /* test polynomials with random Gaussian integer roots */
    for (iter = 0; iter <= 5000; iter++)
    {
        slong num_real_roots;
        int i;
        fmpq_t a, max_root;
        slong bound;
        fmpz_poly_t p,q;

        fmpz_poly_init(p);
        fmpz_poly_init(q);
        fmpq_init(a);
        fmpq_init(max_root);


        fmpz_poly_one(p);

        /* generate random rational roots a/b: (b x - a)*/
        fmpq_set_si(max_root, -1, 1);
        num_real_roots = (slong)n_randint(state, 50);
        for (i = 0; i < num_real_roots; i++)
        {
            fmpq_randtest(a, state, 1 + n_randint(state, 100));
            if(fmpq_cmp(max_root, a) < 0)
                fmpq_set(max_root, a);
            fmpq_neg(a, a);
            fmpz_poly_set_coeff_fmpz(q, 0, fmpq_numref(a));
            fmpz_poly_set_coeff_fmpz(q, 1, fmpq_denref(a));

            fmpz_poly_mul(p, p, q);
        }

        if (fmpz_poly_degree(p) != num_real_roots)
        {
            flint_printf("FAIL:\n");
            flint_abort();
        }

        fmpz_poly_randtest_no_real_root(q, state, 1 + (slong)n_randint(state, 50), 10 + n_randint(state, 100));
        fmpz_poly_mul(p, p, q);

        bound = fmpz_poly_positive_root_upper_bound_2exp(p);

        if (bound == WORD_MIN)
        {
           if (fmpq_sgn(max_root) > 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("got bound = -1 for p = ");
                fmpz_poly_print(p);
                flint_printf("\n");
                flint_printf("but max_root = "); fmpq_print(max_root); flint_printf("\n");
                flint_abort();
            }
        }
        else
        {
            fmpq_one(a);

            if (bound > 0)
                fmpq_mul_2exp(a, a, (ulong)bound);
            else if (bound < 0)
                fmpq_div_2exp(a, a, (ulong)-bound);

            if (fmpq_cmp(a, max_root) < 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("p = "); fmpz_poly_print(p); flint_printf("\n");
                flint_printf("max_real_root = "); fmpq_print(max_root); flint_printf("\n");
                flint_printf("got bound = %wd\n", bound);
                flint_abort();
            }
        }

        fmpq_clear(a);
        fmpq_clear(max_root);
        fmpz_poly_clear(p);
        fmpz_poly_clear(q);
    }


    TEST_FUNCTION_END(state);
}

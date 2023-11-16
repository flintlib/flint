/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "fmpq_poly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_fmpz_poly_composed_op, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, j, k;
        fmpz *r1, *r2;
        fmpz_poly_t A, B, C;
        fmpq_t x, y;
        int op;

        m = 1 + n_randint(state, 8);
        n = 1 + n_randint(state, 8);
        op = n_randint(state, 4);
        i = n_randint(state, m);
        j = n_randint(state, n);

        op = n_randint(state, 4);

        r1 = _fmpz_vec_init(m);
        r2 = _fmpz_vec_init(n);
        fmpz_poly_init(A);
        fmpz_poly_init(B);
        fmpz_poly_init(C);
        fmpq_init(x);
        fmpq_init(y);

        _fmpz_vec_randtest(r1, state, m, 10);
        _fmpz_vec_randtest(r2, state, n, 10);

        /* Don't divide by zero */
        if (op == 3)
        {
            for (k = 0; k < n; k++)
            {
                if (fmpz_is_zero(r2 + k))
                    fmpz_one(r2 + k);
            }
        }

        fmpz_poly_product_roots_fmpz_vec(A, r1, m);
        fmpz_poly_product_roots_fmpz_vec(B, r2, n);

        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(fmpq_numref(x), state, 10);
            fmpz_poly_scalar_mul_fmpz(A, A, fmpq_numref(x));
            fmpz_randtest_not_zero(fmpq_numref(x), state, 10);
            fmpz_poly_scalar_mul_fmpz(B, B, fmpq_numref(x));
        }

        /* output noise */
        fmpz_poly_randtest(C, state, 10, 100);

        qqbar_fmpz_poly_composed_op(C, A, B, op);

        if (op == 0)
        {
            fmpz_add(fmpq_numref(x), r1 + i, r2 + j);
            fmpz_one(fmpq_denref(x));
        }
        else if (op == 1)
        {
            fmpz_sub(fmpq_numref(x), r1 + i, r2 + j);
            fmpz_one(fmpq_denref(x));
        }
        else if (op == 2)
        {
            fmpz_mul(fmpq_numref(x), r1 + i, r2 + j);
            fmpz_one(fmpq_denref(x));
        }
        else
        {
            fmpq_set_fmpz_frac(x, r1 + i, r2 + j);
        }

        fmpz_poly_evaluate_fmpq(y, C, x);

        if (!fmpq_is_zero(y))
        {
            flint_printf("FAIL!\n");
            flint_printf("op = %d\n", op);
            flint_printf("A = "); fmpz_poly_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpz_poly_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpz_poly_print(C); flint_printf("\n\n");
            flint_printf("r1 = "); fmpz_print(r1 + i); flint_printf("\n\n");
            flint_printf("r2 = "); fmpz_print(r2 + j); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_abort();
        }

        _fmpz_vec_clear(r1, m);
        _fmpz_vec_clear(r2, n);
        fmpz_poly_clear(A);
        fmpz_poly_clear(B);
        fmpz_poly_clear(C);
        fmpq_clear(x);
        fmpq_clear(y);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "qqbar.h"

TEST_FUNCTION_START(qqbar_fmpz_poly_symmetric_composed_op, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong d, m, i, j, k;
        fmpz *r, *s;
        fmpz_poly_t A, C, D;
        fmpz_t c;
        int op;

        d = 2 + n_randint(state, 8);
        m = d * (d - 1) / 2;
        op = n_randint(state, 3);

        r = _fmpz_vec_init(d);
        s = _fmpz_vec_init(m);
        fmpz_poly_init(A);
        fmpz_poly_init(C);
        fmpz_poly_init(D);
        fmpz_init(c);

        _fmpz_vec_randtest(r, state, d, 1 + n_randint(state, 10));
        fmpz_poly_product_roots_fmpz_vec(A, r, d);

        if (n_randint(state, 2))
        {
            fmpz_randtest_not_zero(c, state, 10);
            fmpz_poly_scalar_mul_fmpz(A, A, c);
        }

        for (i = k = 0; i < d; i++)
        {
            for (j = i + 1; j < d; j++, k++)
            {
                if (op == 0)
                    fmpz_add(s + k, r + i, r + j);
                else if (op == 1)
                {
                    fmpz_sub(s + k, r + i, r + j);
                    fmpz_mul(s + k, s + k, s + k);
                }
                else
                    fmpz_mul(s + k, r + i, r + j);
            }
        }

        fmpz_poly_product_roots_fmpz_vec(D, s, m);

        /* output noise */
        fmpz_poly_randtest(C, state, 10, 100);

        qqbar_fmpz_poly_symmetric_composed_op(C, A, op);

        fmpz_poly_primitive_part(C, C);
        if (fmpz_poly_length(C) > 0 && fmpz_sgn(fmpz_poly_lead(C)) < 0)
            fmpz_poly_neg(C, C);

        if (!fmpz_poly_equal(C, D))
        {
            flint_printf("FAIL!\n");
            flint_printf("op = %d\n", op);
            flint_printf("A = "); fmpz_poly_print(A); flint_printf("\n\n");
            flint_printf("C = "); fmpz_poly_print(C); flint_printf("\n\n");
            flint_printf("D = "); fmpz_poly_print(D); flint_printf("\n\n");
            flint_abort();
        }

        _fmpz_vec_clear(r, d);
        _fmpz_vec_clear(s, m);
        fmpz_poly_clear(A);
        fmpz_poly_clear(C);
        fmpz_poly_clear(D);
        fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}

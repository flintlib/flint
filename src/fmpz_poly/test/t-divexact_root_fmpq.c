/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(fmpz_poly_divexact_root_fmpq, state)
{
    int i, result;

    /* Compare with standard divrem */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, Q, D, DQ;
        fmpq_t c;
        slong n, b;

        n = n_randint(state, 100);
        b = n_randint(state, 200) + 1;

        fmpq_init(c);
        fmpz_poly_init(P);
        fmpz_poly_init(Q);
        fmpz_poly_init(D);
        fmpz_poly_init(DQ);

        fmpz_poly_randtest(P, state, n, b);
        fmpq_randtest(c, state, b);

        fmpz_poly_set_coeff_fmpz(D, 0, fmpq_numref(c));
        fmpz_poly_neg(D, D);
        fmpz_poly_set_coeff_fmpz(D, 1, fmpq_denref(c));
        fmpz_poly_mul(P, P, D);

        fmpz_poly_divexact_root_fmpq(Q, P, c);

        fmpz_poly_div_basecase(DQ, P, D);

        result = fmpz_poly_equal(Q, DQ);

        if (!result)
        {
            flint_printf("FAIL!\n");
            flint_printf("P:\n"); fmpz_poly_print(P); flint_printf("\n\n");
            flint_printf("Q:\n"); fmpz_poly_print(Q); flint_printf("\n\n");
            flint_printf("D:\n"); fmpz_poly_print(D); flint_printf("\n\n");
            flint_printf("DQ:\n"); fmpz_poly_print(DQ); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(c);
        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);
        fmpz_poly_clear(D);
        fmpz_poly_clear(DQ);
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_poly_t P, D, Q1, Q2;
        fmpq_t c;
        slong n, b;

        n = n_randint(state, 100);
        b = n_randint(state, 200) + 1;

        fmpq_init(c);
        fmpz_poly_init(P);
        fmpz_poly_init(D);
        fmpz_poly_init(Q1);
        fmpz_poly_init(Q2);

        fmpq_randtest(c, state, b);
        fmpz_poly_randtest(P, state, n, b);
        fmpz_poly_set_coeff_fmpz(D, 0, fmpq_numref(c));
        fmpz_poly_neg(D, D);
        fmpz_poly_set_coeff_fmpz(D, 1, fmpq_denref(c));
        fmpz_poly_mul(P, P, D);

        fmpz_poly_set(Q2, P);

        fmpz_poly_divexact_root_fmpq(Q1, P, c);
        fmpz_poly_divexact_root_fmpq(Q2, Q2, c);

        result = fmpz_poly_equal(Q1, Q2);

        if (!result)
        {
            flint_printf("FAIL (aliasing)!\n");
            flint_printf("P:\n"); fmpz_poly_print(P); flint_printf("\n\n");
            flint_printf("Q1:\n"); fmpz_poly_print(Q1); flint_printf("\n\n");
            flint_printf("Q2:\n"); fmpz_poly_print(Q2); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(c);
        fmpz_poly_clear(P);
        fmpz_poly_clear(D);
        fmpz_poly_clear(Q1);
        fmpz_poly_clear(Q2);
    }

    TEST_FUNCTION_END(state);
}

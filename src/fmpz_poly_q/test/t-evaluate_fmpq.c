/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_poly_q.h"

TEST_FUNCTION_START(fmpz_poly_q_evaluate_fmpq, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        int ans1, ans2;
        fmpq_t a, b;
        fmpz_t num, den;
        fmpz_poly_q_t f;

        fmpq_init(a);
        fmpq_init(b);
        fmpz_init(num);
        fmpz_init(den);
        fmpz_poly_q_init(f);
        fmpz_poly_q_randtest(f, state, n_randint(state, 10), 10, n_randint(state, 10), 10);

        fmpz_randtest(num, state, 50);
        fmpz_randtest_not_zero(den, state, 50);
        fmpz_set(fmpq_numref(a), num);
        fmpz_set(fmpq_denref(a), den);
        fmpq_canonicalise(a);

        ans1 = fmpz_poly_q_evaluate_fmpq(b, f, a);
        ans2 = fmpz_poly_q_evaluate_fmpq(a, f, a);

        result = (ans1 == ans2) && fmpq_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "), fmpz_poly_q_print(f), flint_printf("\n");
            flint_printf("num = "), fmpz_print(num), flint_printf("\n");
            flint_printf("den = "), fmpz_print(den), flint_printf("\n");
            flint_printf("a = "); fmpq_print(a); flint_printf("\n");
            flint_printf("b = "); fmpq_print(b); flint_printf("\n");
            flint_printf("ans1 = %d\n", ans1);
            flint_printf("ans2 = %d\n", ans2);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(a);
        fmpq_clear(b);
        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_poly_q_clear(f);
    }

    TEST_FUNCTION_END(state);
}

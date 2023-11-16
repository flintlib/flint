/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_exp_sum_bs_powtab, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong N;
        fmpz_t x, T, Q;
        fmpq_t S, V;
        flint_bitcnt_t Qexp, r;

        fmpz_init(x);
        fmpz_init(T);
        fmpz_init(Q);
        fmpq_init(S);
        fmpq_init(V);

        if (n_randint(state, 100) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));
        else
            flint_set_num_threads(1);

        N = 1 + n_randint(state, 300);
        r = n_randint(state, 10);
        fmpz_randtest(x, state, 80);

        _arb_exp_sum_bs_simple(T, Q, &Qexp, x, r, N);
        fmpq_set_fmpz_frac(S, T, Q);
        fmpq_div_2exp(S, S, Qexp);

        _arb_exp_sum_bs_powtab(T, Q, &Qexp, x, r, N);
        fmpq_set_fmpz_frac(V, T, Q);
        fmpq_div_2exp(V, V, Qexp);

        if (!fmpq_equal(S, V))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wd\n\n", N);
            flint_printf("r = %wu\n\n", r);
            flint_printf("x = "); fmpz_print(x); flint_printf("\n\n");
            flint_printf("T = "); fmpz_print(T); flint_printf("\n\n");
            flint_printf("Q = "); fmpz_print(T); flint_printf("\n\n");
            flint_printf("V = "); fmpq_print(V); flint_printf("\n\n");
            flint_printf("S = "); fmpq_print(S); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(x);
        fmpz_clear(T);
        fmpz_clear(Q);
        fmpq_clear(S);
        fmpq_clear(V);
    }

    TEST_FUNCTION_END(state);
}

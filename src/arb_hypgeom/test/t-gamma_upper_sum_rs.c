/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_gamma_upper_sum_rs, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        ulong p, q;
        arb_t z, r1, r2, s, t, u;
        slong k, N, prec;

        p = n_randint(state, 1000);
        q = 1 + n_randint(state, 1000);

        N = n_randint(state, 100);
        prec = 2 + n_randint(state, 500);

        fmpq_init(a);
        fmpq_set_si(a, p, q);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);
        arb_init(s);
        arb_init(t);
        arb_init(u);

        arb_randtest(z, state, prec, 10);

        _arb_hypgeom_gamma_upper_sum_rs_1(r1, p, q, z, N, prec);

        arb_zero(s);
        arb_one(t);
        arb_set_fmpq(u, a, prec);
        for (k = 0; k < N; k++)
        {
            arb_add(s, s, t, prec);
            arb_div(t, t, z, prec);
            arb_mul(t, t, u, prec);
            arb_add_ui(u, u, 1, prec);
        }
        arb_set(r2, s);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("N = %wd\n\n", N);
            flint_printf("a = "); fmpq_print(a); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 100, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_clear(a);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);
        arb_clear(s);
        arb_clear(t);
        arb_clear(u);
    }

    TEST_FUNCTION_END(state);
}

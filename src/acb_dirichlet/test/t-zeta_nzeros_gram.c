/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_zeta_nzeros_gram, state)
{
    slong iter;

    for (iter = 0; iter < 130 + 20 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t t, x;
        fmpz_t N, n;
        slong prec1, prec2;

        arb_init(t);
        arb_init(x);
        fmpz_init(n);
        fmpz_init(N);

        if (iter < 130)
        {
            fmpz_set_si(n, iter - 1);
        }
        else
        {
            fmpz_randtest_unsigned(n, state, 20);
            fmpz_add_ui(n, n, 129);
        }

        prec1 = 2 + n_randtest(state) % 100;
        prec2 = 2 + n_randtest(state) % 100;

        acb_dirichlet_zeta_nzeros_gram(N, n);

        acb_dirichlet_gram_point(t, n, NULL, NULL, prec1);
        acb_dirichlet_zeta_nzeros(x, t, prec2);

        if (!arb_contains_fmpz(x, N))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("n = "); fmpz_print(n);
            flint_printf("   prec1 = %wd  prec2 = %wd\n\n", prec1, prec2);
            flint_printf("N = "); fmpz_print(N); flint_printf("\n\n");
            flint_printf("t = "); arb_printn(t, 100, 0); flint_printf("\n\n");
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(t);
        arb_clear(x);
        fmpz_clear(n);
        fmpz_clear(N);
    }

    TEST_FUNCTION_END(state);
}

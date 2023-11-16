/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

void acb_dirichlet_stieltjes_integral(acb_t res, const fmpz_t n, const acb_t a, slong prec);
void acb_dirichlet_stieltjes_em(acb_t res, const fmpz_t n, const acb_t a, slong prec);

TEST_FUNCTION_START(acb_dirichlet_stieltjes, state)
{
    slong iter;

    for (iter = 0; iter < 250 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t b1, b2, a;
        fmpz_t n;
        slong prec1, prec2, acc1, acc2;
        int alg1, alg2;

        fmpz_init(n);
        acb_init(b1);
        acb_init(b2);
        acb_init(a);

        if (n_randint(state, 10) == 0)
        {
            acb_one(a);
            fmpz_randtest(n, state, 1 + n_randint(state, 100));
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 100);
            prec2 = 2 + n_randint(state, 100);
        }
        else if (n_randint(state, 7) != 0)
        {
            acb_one(a);
            fmpz_randtest(n, state, 2 + n_randint(state, 12));
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 400);
            prec2 = 2 + n_randint(state, 400);
        }
        else
        {
            acb_randtest_precise(a, state, 400, 3);
            fmpz_randtest(n, state, 7);
            fmpz_abs(n, n);
            prec1 = 2 + n_randint(state, 400);
            prec2 = 2 + n_randint(state, 400);
        }

        alg1 = n_randint(state, 3);
        alg2 = n_randint(state, 3);

        if (fmpz_cmp_ui(n, 40) >= 0)
        {
            alg1 = 1 + n_randint(state, 2);
            alg2 = 1 + n_randint(state, 2);
        }
        else
        {
            alg1 = n_randint(state, 3);
            alg2 = n_randint(state, 3);
        }

        if (alg1 == 0)
            acb_dirichlet_stieltjes_em(b1, n, a, prec1);
        else if (alg1 == 1)
            acb_dirichlet_stieltjes_integral(b1, n, a, prec1);
        else
            acb_dirichlet_stieltjes(b1, n, a, prec1);

        if (alg2 == 0)
            acb_dirichlet_stieltjes_em(b2, n, a, prec2);
        else if (alg2 == 1)
            acb_dirichlet_stieltjes_integral(b2, n, a, prec2);
        else
            acb_dirichlet_stieltjes(b2, n, a, prec2);

        if (!acb_overlaps(b1, b2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd, alg1 = %d, alg2 = %d\n\n", iter, alg1, alg2);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("a = "); acb_printn(a, 100, 0); flint_printf("\n\n");
            flint_printf("b1 = "); acb_printn(b1, 1000, 0); flint_printf("\n\n");
            flint_printf("b2 = "); acb_printn(b2, 1000, 0); flint_printf("\n\n");
            flint_abort();
        }

        acc1 = acb_rel_accuracy_bits(b1);
        acc2 = acb_rel_accuracy_bits(b2);

        if (acb_is_one(a) && (acc1 < prec1 - 10 || acc2 < prec2 - 10))
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec1 = %wd, acc1 = %wd\n", prec1, acc1);
            flint_printf("prec2 = %wd, acc2 = %wd\n", prec2, acc2);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n\n");
            flint_printf("b1 = "); acb_printn(b1, 500, 0); flint_printf("\n\n");
            flint_printf("b2 = "); acb_printn(b2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(b1);
        acb_clear(b2);
        acb_clear(a);
        fmpz_clear(n);
    }

    TEST_FUNCTION_END(state);
}

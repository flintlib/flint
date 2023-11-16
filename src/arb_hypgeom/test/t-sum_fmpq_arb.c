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
#include "fmpq_vec.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_sum_fmpq_arb, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t s1, s2, s3, z;
        fmpq *a, *b;
        slong alen, blen, N, prec;
        int reciprocal;

        alen = n_randint(state, 5);
        blen = n_randint(state, 5);
        N = n_randint(state, 100);

        arb_init(s1);
        arb_init(s2);
        arb_init(s3);
        arb_init(z);
        a = _fmpq_vec_init(alen);
        b = _fmpq_vec_init(blen);

        prec = 2 + n_randint(state, 500);
        reciprocal = n_randint(state, 2);

        if (n_randint(state, 10) == 0)
            arb_randtest_special(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        else
            arb_randtest(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        arb_randtest_special(s1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest_special(s2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        _fmpq_vec_randtest(a, state, alen, 1 + n_randint(state, 100));
        _fmpq_vec_randtest(b, state, blen, 1 + n_randint(state, 100));

        arb_hypgeom_sum_fmpq_arb_forward(s1, a, alen, b, blen, z, reciprocal, N, prec);
        arb_hypgeom_sum_fmpq_arb_rs(s2, a, alen, b, blen, z, reciprocal, N, prec);
        arb_hypgeom_sum_fmpq_arb_bs(s3, a, alen, b, blen, z, reciprocal, N, prec);

        if (!arb_overlaps(s1, s2) || !arb_overlaps(s1, s3))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("N = %wd\n\n", N);
            flint_printf("a = "); _fmpq_vec_print(a, alen); flint_printf("\n\n");
            flint_printf("b = "); _fmpq_vec_print(b, blen); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 100, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 100, 0); flint_printf("\n\n");
            flint_printf("s3 = "); arb_printn(s3, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(s1);
        arb_clear(s2);
        arb_clear(s3);
        arb_clear(z);
        _fmpq_vec_clear(a, alen);
        _fmpq_vec_clear(b, blen);
    }

    TEST_FUNCTION_END(state);
}

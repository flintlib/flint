/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"

TEST_FUNCTION_START(fmpq_mat_trace, state)
{
    slong i;

    /* Test trace(AB) = trace(BA) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mat_t A, B, AB, BA;
        fmpq_t trab, trba;
        slong m, n;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, m);
        fmpq_mat_init(AB, m, m);
        fmpq_mat_init(BA, n, n);

        fmpq_init(trab);
        fmpq_init(trba);

        fmpq_mat_randtest(A, state, 1 + n_randint(state, 100));
        fmpq_mat_randtest(B, state, 1 + n_randint(state, 100));

        fmpq_mat_mul(AB, A, B);
        fmpq_mat_mul(BA, B, A);

        fmpq_mat_trace(trab, AB);
        fmpq_mat_trace(trba, BA);

        if (!fmpq_equal(trab, trba))
        {
            flint_printf("FAIL:\n");
            fmpq_mat_print(A), flint_printf("\n");
            fmpq_mat_print(B), flint_printf("\n");
            fmpq_mat_print(AB), flint_printf("\n");
            fmpq_mat_print(BA), flint_printf("\n");
            flint_printf("tr(AB): "),  fmpq_print(trab),    flint_printf("\n");
            flint_printf("tr(BA): "),  fmpq_print(trba),    flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(AB);
        fmpq_mat_clear(BA);
        fmpq_clear(trab);
        fmpq_clear(trba);
    }

    TEST_FUNCTION_END(state);
}

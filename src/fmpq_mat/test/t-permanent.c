/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "gr.h"
#include "gr_mat.h"

TEST_FUNCTION_START(fmpq_mat_permanent, state)
{
    slong ix;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpq_mat_t A;
        fmpq_t a, b;
        slong n, bits;
        gr_ctx_t ctx;

        gr_ctx_init_fmpq(ctx);

        n = n_randint(state, 4);
        bits = 1 + n_randint(state, 100);

        fmpq_init(a);
        fmpq_init(b);

        fmpq_mat_init(A, n, n);
        fmpq_mat_randtest(A, state, bits);

        if (!fmpq_mat_permanent(a, A))
        {
            flint_printf("FAIL: unsuccessful with n = %wd\n", n);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_mat_permanent_generic(b, (gr_mat_struct *) A, ctx));

        if (!fmpq_equal(a, b))
        {
            flint_printf("FAIL:\n");
            fmpq_mat_print(A), flint_printf("\n");
            flint_printf("expected: "),  fmpq_print(b),    flint_printf("\n");
            flint_printf("computed: "), fmpq_print(a), flint_printf("\n");
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_clear(a);
        fmpq_clear(b);
    }

    TEST_FUNCTION_END(state);
}

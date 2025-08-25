/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

TEST_FUNCTION_START(fmpz_mat_permanent, state)
{
    slong ix;

    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        fmpz_mat_t A;
        fmpz_t a, b;
        slong n, bits;
        gr_ctx_t ctx;

        gr_ctx_init_fmpz(ctx);

        flint_set_num_threads(1 + n_randint(state, 3));

        if (n_randint(state, 100) == 0)
        {
            n = n_randint(state, 16);
            bits = 1 + n_randint(state, 64);
        }
        else
        {
            n = n_randint(state, 10);
            bits = 1 + n_randint(state, 100);
        }

        fmpz_init(a);
        fmpz_init(b);

        fmpz_mat_init(A, n, n);
        fmpz_mat_randtest(A, state, bits);

        if (!fmpz_mat_permanent(a, A))
        {
            flint_printf("FAIL: unsuccessful with n = %wd\n", n);
            flint_abort();
        }

        GR_MUST_SUCCEED(gr_mat_permanent_generic(b, (gr_mat_struct *) A, ctx));

        if (!fmpz_equal(a, b))
        {
            flint_printf("FAIL:\n");
            fmpz_mat_print_pretty(A), flint_printf("\n");
            flint_printf("expected: "),  fmpz_print(b),    flint_printf("\n");
            flint_printf("computed: "), fmpz_print(a), flint_printf("\n");
            flint_abort();
        }

        fmpz_mat_clear(A);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}

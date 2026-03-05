/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_complex, state)
{
    gr_ctx_t R, C;
    int reps = 10 * flint_test_multiplier();
    int i, test_flags;

    for (i = 0; i < 4; i++)
    {
        test_flags = GR_TEST_ALWAYS_ABLE;
        test_flags = 0;

        switch (i)
        {
            case 0:  gr_ctx_init_fmpz(R); break;
            case 1:  gr_ctx_init_fmpq(R); break;
            case 2:  gr_ctx_init_fmpz_poly(R); break;
            case 3:  gr_ctx_init_real_arb(R, 53); break;
        }

        gr_ctx_init_gr_complex(C, R);
        gr_test_ring(C, reps, test_flags);

        gr_ctx_clear(C);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}

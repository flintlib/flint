/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_fraction, state)
{
    gr_ctx_t R, R2, K;
    int reps = 10 * flint_test_multiplier();
    int i, test_flags, haveR2;

    for (i = 0; i < 7; i++)
    {
        haveR2 = 0;
        test_flags = GR_TEST_ALWAYS_ABLE;

        switch (i)
        {
            case 0:  gr_ctx_init_fmpz(R); break;
            case 1:  gr_ctx_init_fmpq(R); break;
            case 2:  gr_ctx_init_fmpz_poly(R); break;
            case 3:  gr_ctx_init_fmpq_poly(R); break;
            case 4:  gr_ctx_init_real_arb(R, 53); test_flags &= ~GR_TEST_ALWAYS_ABLE; break;
            case 5:  gr_ctx_init_fmpzi(R); break;
            case 6:  gr_ctx_init_fmpzi(R2); haveR2 = 1; gr_ctx_init_gr_poly(R, R2); break;
        }

        gr_ctx_init_gr_fraction(K, R, 0);
        gr_test_ring(K, reps, test_flags);
        gr_ctx_clear(K);

        gr_ctx_init_gr_fraction(K, R, GR_FRACTION_NO_REDUCTION);
        gr_test_ring(K, reps, test_flags);
        gr_ctx_clear(K);

        gr_ctx_init_gr_fraction(K, R, GR_FRACTION_STRONGLY_CANONICAL);
        gr_test_ring(K, reps, test_flags);
        gr_ctx_clear(K);

        gr_ctx_clear(R);
        if (haveR2)
            gr_ctx_clear(R2);
    }

    TEST_FUNCTION_END(state);
}

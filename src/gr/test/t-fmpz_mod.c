/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "gr.h"

TEST_FUNCTION_START(gr_fmpz_mod, state)
{
    gr_ctx_t ZZn;
    fmpz_t n;
    slong iter;
    int flags = GR_TEST_ALWAYS_ABLE;

    fmpz_init(n);

    for (iter = 0; iter < 100; iter++)
    {
        fmpz_randtest_not_zero(n, state, 200);
        fmpz_abs(n, n);
        gr_ctx_init_fmpz_mod(ZZn, n);
        gr_test_ring(ZZn, 100, flags);
        gr_ctx_clear(ZZn);
    }

    /* test huge preinvn code */
    {
        fmpz_randbits(n, state, 72000);
        fmpz_abs(n, n);
        gr_ctx_init_fmpz_mod(ZZn, n);
        gr_test_ring(ZZn, 10, flags);
        gr_ctx_clear(ZZn);
    }

    fmpz_clear(n);

    TEST_FUNCTION_END(state);
}

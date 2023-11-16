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

TEST_FUNCTION_START(gr_fmpz, state)
{
    gr_ctx_t ZZ;
    int flags = GR_TEST_ALWAYS_ABLE;

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 1000;
    gr_test_ring(ZZ, 100000, flags);

    {
        fmpz * a, *b;
        a = gr_heap_init(ZZ);
        b = gr_heap_init(ZZ);

        fmpz_set_str(a, "1000000000000000000000", 10);
        fmpz_set_str(b, "1000000000000000000001", 10);

        if (gr_sub(b, b, a, ZZ) != GR_SUCCESS || gr_is_one(b, ZZ) != T_TRUE)
            flint_abort();

        gr_heap_clear(a, ZZ);
        gr_heap_clear(b, ZZ);
    }

    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}

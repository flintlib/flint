/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "ca.h"

TEST_FUNCTION_START(ca_get_str, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x;
        char * s;
        slong slen;

        ca_ctx_init(ctx);
        ca_init(x, ctx);

        /* just check that it doesn't crash */
        ca_randtest_special(x, state, 10, 100, ctx);
        s = ca_get_str(x, ctx);
        slen = strlen(s);
        (void)slen;
        flint_free(s);

        ca_clear(x, ctx);
        ca_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

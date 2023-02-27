/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "ca.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_str...");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
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
        slen = slen;
        flint_free(s);

        ca_clear(x, ctx);
        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


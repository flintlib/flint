/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("get_fexpr....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_t x;
        fexpr_t f;

        ca_ctx_init(ctx);
        ca_init(x, ctx);
        fexpr_init(f);

        /* Just check that this runs, for now... */
        ca_randtest_special(x, state, 5, 5, ctx);
        ca_get_fexpr(f, x, 0, ctx);

        fexpr_clear(f);
        ca_clear(x, ctx);
        ca_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("cache_insert....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * calcium_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_ext_t ext;
        ca_ext_ptr ext2;
        ca_ext_cache_t cache;
        qqbar_t x;
        slong i, len;

        ca_ctx_init(ctx);
        qqbar_init(x);

        len = n_randint(state, 1000);

        ca_ext_cache_init(cache, ctx);

        for (i = 0; i < len; i++)
        {
            qqbar_set_ui(x, n_randint(state, 100));
            ca_ext_init_qqbar(ext, x, ctx);

            ext2 = ca_ext_cache_insert(cache, ext, ctx);

            if (!ca_ext_equal_repr(ext2, ext, ctx))
            {
                flint_printf("FAIL\n\n");
                flint_printf("x = "); qqbar_print(x); flint_printf("\n\n");
                flint_printf("ext = "); ca_ext_print(ext, ctx); flint_printf("\n\n");
                flint_printf("ext2 = "); ca_ext_print(ext2, ctx); flint_printf("\n\n");
                flint_abort();
            }

            ca_ext_clear(ext, ctx);
        }

        ca_ext_cache_clear(cache, ctx);
        ca_ctx_clear(ctx);
        qqbar_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


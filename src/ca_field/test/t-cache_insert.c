/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

static int _ca_field_equal_ext(const ca_field_t K, ca_ext_struct ** x, slong len, ca_ctx_t ctx)
{
    slong i;

    if (len != CA_FIELD_LENGTH(K))
        return 0;

    for (i = 0; i < len; i++)
        if (CA_FIELD_EXT_ELEM(K, i) != x[i])
            return 0;

    return 1;
}

TEST_FUNCTION_START(ca_field_cache_insert, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ca_ctx_t ctx;
        ca_ext_struct ** ext;
        ca_ext_struct * all_ext;
        ca_field_cache_t cache;
        ca_field_struct *K, *K2;
        qqbar_t x;
        ca_t v;
        slong i, j, steps, len, ext_len;

        ca_ctx_init(ctx);
        qqbar_init(x);
        ca_init(v, ctx);

        ext_len = 1 + n_randint(state, 10);
        steps = n_randint(state, 100);

        all_ext = flint_malloc(sizeof(ca_ext_struct) * ext_len);
        for (j = 0; j < ext_len; j++)
        {
            if (n_randint(state, 10) != 0)
            {
                ca_set_ui(v, j + 2, ctx);
                ca_ext_init_fx(all_ext + j, CA_RiemannZeta, v, ctx);
            }
            else
            {
                qqbar_set_ui(x, j);
                ca_ext_init_qqbar(all_ext + j, x, ctx);
            }
        }

        ca_field_cache_init(cache, ctx);

        for (i = 0; i < steps; i++)
        {
            len = n_randint(state, 10);

            ext = flint_malloc(len * sizeof(ca_ext_struct *));

            for (j = 0; j < len; j++)
                ext[j] = all_ext + n_randint(state, ext_len);

            K = ca_field_cache_insert_ext(cache, ext, len, ctx);
            K2 = ca_field_cache_insert_ext(cache, ext, len, ctx);

            if (!_ca_field_equal_ext(K, ext, len, ctx) || K != K2)
            {
                flint_printf("FAIL\n\n");
                flint_printf("K = "); ca_field_print(K, ctx); flint_printf("\n\n");
                flint_printf("ext = \n\n");
                for (j = 0; j < len; j++)
                {
                    ca_ext_print(ext[j], ctx); flint_printf("\n\n");
                }
                flint_abort();
            }

            flint_free(ext);
        }

        for (j = 0; j < ext_len; j++)
            ca_ext_clear(all_ext + j, ctx);

        flint_free(all_ext);

        ca_field_cache_clear(cache, ctx);
        ca_clear(v, ctx);
        ca_ctx_clear(ctx);
        qqbar_clear(x);
    }

    TEST_FUNCTION_END(state);
}

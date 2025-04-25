/*
    Copyright (C) 2025 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_vec.h"

TEST_FUNCTION_START(gr_vec_sort, state)
{
    for (slong iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_vec_t vec0, vec1, vec2;

        gr_ctx_init_random(ctx, state);
        slong len = n_randint(state, 30);
        gr_vec_init(vec0, len, ctx);
        gr_vec_init(vec1, len, ctx);
        gr_vec_init(vec2, 0, ctx);

        int status = GR_SUCCESS;

        status |= _gr_vec_randtest(vec0->entries, state, len, ctx);

        status |= gr_vec_set(vec1, vec0, ctx);
        _gr_vec_shuffle(vec1->entries, state, len, ctx);

        status |= gr_vec_sort(vec2, vec0, ctx);
        status |= gr_vec_sort(vec1, vec1, ctx);

        if (status == GR_SUCCESS)
        {
            if(_gr_vec_equal(vec1->entries, vec2->entries, len, ctx) == T_FALSE)
                status = GR_TEST_FAIL;

            for (slong i = 0; i < len - 1; i++)
                if (gr_le(gr_vec_entry_ptr(vec1, i, ctx),
                          gr_vec_entry_ptr(vec1, i + 1, ctx), ctx) == T_FALSE)
                    status = GR_TEST_FAIL;
        }

        for (slong i = 0; i < len - 1; i++)
            if (_gr_vec_contains(vec0->entries, len,
                                 gr_vec_entry_ptr(vec1, i, ctx),
                                 ctx) == T_FALSE)
                status = GR_TEST_FAIL;

        if (status == GR_TEST_FAIL)
        {
            flint_printf("FAIL\n");
            printf("vec0 = "); _gr_vec_print(vec0, len, ctx); printf("\n");
            flint_abort();
        }

        gr_vec_clear(vec2, ctx);
        gr_vec_clear(vec1, ctx);
        gr_vec_clear(vec0, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}


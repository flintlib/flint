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
#include "perm.h"

TEST_FUNCTION_START(gr_vec_permute, state)
{
    for (slong iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_vec_t vec, out;

        gr_ctx_init_random(ctx, state);
        slong len = n_randint(state, 30);
        gr_vec_init(vec, len, ctx);
        gr_vec_init(out, 0, ctx);

        slong * perm = _perm_init(len);
        _perm_randtest(perm, len, state);

        int status = GR_SUCCESS;

        status |= _gr_vec_randtest(vec->entries, state, len, ctx);

        status |= gr_vec_permute(out, vec, perm, ctx);

        if (status == GR_SUCCESS)
        {
            for (slong k = 0; k < len; k++)
                if (gr_equal(gr_vec_entry_ptr(vec, k, ctx),
                             gr_vec_entry_ptr(out, perm[k], ctx),
                             ctx) == T_FALSE)
                    status = GR_TEST_FAIL;
        }

        _perm_inv(perm, perm, len);
        status |= gr_vec_permute(out, out, perm, ctx);

        if (status == GR_SUCCESS &&
            _gr_vec_equal(out->entries, vec->entries, len, ctx) == T_FALSE)
            status = GR_TEST_FAIL;

        if (status == GR_TEST_FAIL)
        {
            flint_printf("FAIL\n");
            printf("vec = "); _gr_vec_print(vec->entries, len, ctx); printf("\n");
            flint_abort();
        }

        _perm_clear(perm);
        gr_vec_clear(out, ctx);
        gr_vec_clear(vec, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

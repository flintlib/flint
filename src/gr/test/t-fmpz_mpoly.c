/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "mpoly.h"
#include "gr.h"
#include "gr_vec.h"

TEST_FUNCTION_START(gr_fmpz_mpoly, state)
{
    gr_ctx_t ZZxy;
    slong iter;
    int flags = 0;

    for (iter = 0; iter < 10; iter++)
    {
        gr_ctx_init_fmpz_mpoly(ZZxy, n_randint(state, 3), mpoly_ordering_randtest(state));
        ZZxy->size_limit = 100;
        gr_test_ring(ZZxy, 100, flags);
        gr_ctx_clear(ZZxy);
    }

    /* generator names can be set repeatedly */
    {
        const char * names1[] = {"x", "y"};
        const char * names2[] = {"alpha", "beta_1"};
        gr_vec_t gens;
        char * str;

        gr_ctx_init_fmpz_mpoly(ZZxy, 2, ORD_LEX);
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(ZZxy, names1));
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(ZZxy, names2));
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(ZZxy, names1));
        GR_MUST_SUCCEED(gr_ctx_set_gen_names(ZZxy, names2));
        gr_vec_init(gens, 0, ZZxy);
        GR_MUST_SUCCEED(gr_gens(gens, ZZxy));
        GR_MUST_SUCCEED(gr_get_str(&str, gr_vec_entry_srcptr(gens, 1, ZZxy), ZZxy));
        if (strcmp(str, "beta_1") != 0)
        {
            flint_printf("FAIL: gen names after repeated gr_ctx_set_gen_names: %s\n", str);
            flint_abort();
        }
        flint_free(str);
        gr_vec_clear(gens, ZZxy);
        gr_ctx_clear(ZZxy);
    }

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "gr_tower.h"

TEST_FUNCTION_START(gr_tower_flat, state)
{
    gr_ctx_t QQ;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_tower_t T;
        gr_ctx_t K, F;
        qqbar_t x;
        slong k, n = 1 + n_randint(state, 3);
        gr_ptr a, b, c, fa, fb, fc;

        qqbar_init(x);
        gr_tower_init(T, QQ);

        for (k = 0; k < n; k++)
        {
            if (n_randint(state, 2))
            {
                qqbar_randtest(x, state, 1 + n_randint(state, 3), 6);
                GR_MUST_SUCCEED(gr_tower_adjoin_qqbar(T, x, NULL));
            }
            else
            {
                gr_ctx_struct * top = gr_tower_field(T);
                gr_ptr e;
                GR_TMP_INIT(e, top);
                GR_MUST_SUCCEED(gr_randtest(e, state, top));
                GR_IGNORE(gr_tower_adjoin_root_ui(T, e, 2 + n_randint(state, 2), NULL));
                GR_TMP_CLEAR(e, top);
            }
        }

        gr_ctx_init_tower_field(K, T);
        gr_ctx_init_tower_field_flat(F, T);

        gr_test_ring(F, 3 * flint_test_multiplier(), GR_TEST_ALWAYS_ABLE);

        /* nested <-> flat round trips and consistency of arithmetic */
        GR_TMP_INIT3(a, b, c, K);
        GR_TMP_INIT3(fa, fb, fc, F);

        GR_MUST_SUCCEED(gr_randtest(a, state, K));
        GR_MUST_SUCCEED(gr_randtest(b, state, K));
        GR_MUST_SUCCEED(gr_set_other(fa, a, K, F));
        GR_MUST_SUCCEED(gr_set_other(fb, b, K, F));

        /* c = a*b + a - 3, in both */
        GR_MUST_SUCCEED(gr_mul(c, a, b, K));
        GR_MUST_SUCCEED(gr_add(c, c, a, K));
        GR_MUST_SUCCEED(gr_sub_ui(c, c, 3, K));
        GR_MUST_SUCCEED(gr_mul(fc, fa, fb, F));
        GR_MUST_SUCCEED(gr_add(fc, fc, fa, F));
        GR_MUST_SUCCEED(gr_sub_ui(fc, fc, 3, F));

        {
            gr_ptr back;
            GR_TMP_INIT(back, K);
            GR_MUST_SUCCEED(gr_tower_flat_get_nested(back, fc, F));
            if (gr_equal(back, c, K) != T_TRUE)
            {
                flint_printf("FAIL: flat vs nested arithmetic\n");
                gr_println(c, K); gr_println(fc, F); gr_println(back, K);
                flint_abort();
            }
            GR_TMP_CLEAR(back, K);
        }

        /* division in flat, compared after conversion */
        if (gr_is_zero(c, K) == T_FALSE)
        {
            gr_ptr back;
            GR_TMP_INIT(back, K);
            GR_MUST_SUCCEED(gr_div(fc, fa, fc, F));
            GR_MUST_SUCCEED(gr_div(c, a, c, K));
            GR_MUST_SUCCEED(gr_tower_flat_get_nested(back, fc, F));
            if (gr_equal(back, c, K) != T_TRUE)
            {
                flint_printf("FAIL: flat division\n");
                flint_abort();
            }
            /* and equality test directly in the flat representation */
            {
                gr_ptr fc2;
                GR_TMP_INIT(fc2, F);
                GR_MUST_SUCCEED(gr_set_other(fc2, c, K, F));
                if (gr_equal(fc, fc2, F) != T_TRUE)
                {
                    flint_printf("FAIL: flat equality with different denominators\n");
                    gr_println(fc, F); gr_println(fc2, F);
                    flint_abort();
                }
                GR_TMP_CLEAR(fc2, F);
            }
            GR_TMP_CLEAR(back, K);
        }

        GR_TMP_CLEAR3(a, b, c, K);
        GR_TMP_CLEAR3(fa, fb, fc, F);
        gr_ctx_clear(F);
        gr_ctx_clear(K);
        gr_tower_clear(T);
        qqbar_clear(x);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}

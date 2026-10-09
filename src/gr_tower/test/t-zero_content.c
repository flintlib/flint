/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_q.h"
#include "test_helpers.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_tower.h"

/*
    The exact zero test removes the monomial content M of a numerator
    (testing y for x = M y) only when M is proved nonzero by its
    enclosure. Here x = s (e^2 - exp(2)) and x = s (e^2 - exp(2) +
    10^-30) in Q(pi, e, exp(2), s) with exp(2) an independent generator
    (the relation is found by the zero test) and s = sqrt(u) for u = 3 +
    pi, u = pi - (a 62-digit approximation) (a tiny positive number), and
    u = exp(2) - e^2 (zero, but not structurally: the root is not
    adjoined, since u is not proved nonzero).
*/

TEST_FUNCTION_START(gr_tower_zero_content, state)
{
    gr_ctx_t QQ;
    int which;

    gr_ctx_init_fmpq(QQ);

    for (which = 0; which < 3; which++)
    {
        gr_tower_t T;
        gr_ctx_t F;
        gr_ctx_struct * top;
        gr_vec_t g, G;
        gr_ptr u, x, y;
        slong j;
        int status;
        truth_t t;

        gr_tower_init(T, QQ);
        GR_MUST_SUCCEED(gr_tower_adjoin_pi(T, NULL));
        for (j = 1; j <= 2; j++)
        {
            fmpz_mpoly_q_t c;
            fmpz_mpoly_ctx_struct * mctx = gr_tower_flat(T)->mctx;
            fmpz_mpoly_q_init(c, mctx);
            fmpz_mpoly_q_set_si(c, j, mctx);
            GR_MUST_SUCCEED(gr_tower_adjoin_exp_flat(T, c, mctx, NULL));
            fmpz_mpoly_q_clear(c, mctx);
        }

        top = gr_tower_field(T);
        gr_vec_init(g, 0, top);
        GR_MUST_SUCCEED(gr_gens(g, top));
        u = gr_heap_init(top);
        if (which == 0)
        {
            GR_MUST_SUCCEED(gr_set_ui(u, 3, top));
            GR_MUST_SUCCEED(gr_add(u, u, gr_vec_entry_ptr(g, 0, top), top));
        }
        else if (which == 1)
        {
            GR_MUST_SUCCEED(gr_set_str(u, "314159265358979323846264338327950288419716939937510582097494459/10^62", top));
            GR_MUST_SUCCEED(gr_sub(u, gr_vec_entry_ptr(g, 0, top), u, top));
        }
        else
        {
            GR_MUST_SUCCEED(gr_sqr(u, gr_vec_entry_ptr(g, 1, top), top));
            GR_MUST_SUCCEED(gr_sub(u, gr_vec_entry_ptr(g, 2, top), u, top));
        }

        status = gr_tower_adjoin_root_ui(T, u, 2, "s");
        gr_heap_clear(u, top);
        gr_vec_clear(g, top);

        if (which == 2)
        {
            if (status == GR_SUCCESS)
            {
                flint_printf("FAIL: the square root of a zero was adjoined\n");
                flint_abort();
            }
            gr_tower_clear(T);
            continue;
        }
        GR_MUST_SUCCEED(status);

        gr_ctx_init_tower_field_flat(F, T);
        gr_vec_init(G, 0, F);
        GR_MUST_SUCCEED(gr_gens(G, F));
        x = gr_heap_init(F);
        y = gr_heap_init(F);

        /* s (e^2 - exp(2)) = 0 */
        GR_MUST_SUCCEED(gr_sqr(x, gr_vec_entry_ptr(G, 1, F), F));
        GR_MUST_SUCCEED(gr_sub(x, x, gr_vec_entry_ptr(G, 2, F), F));
        GR_MUST_SUCCEED(gr_mul(y, x, gr_vec_entry_ptr(G, 3, F), F));
        t = gr_is_zero(y, F);
        if (t != T_TRUE)
        {
            flint_printf("FAIL: zero (%d), case %d\n", (int) t, which);
            flint_abort();
        }

        /* s (e^2 - exp(2) + 10^-30) != 0 */
        GR_MUST_SUCCEED(gr_set_str(y, "1/10^30", F));
        GR_MUST_SUCCEED(gr_add(y, x, y, F));
        GR_MUST_SUCCEED(gr_mul(y, y, gr_vec_entry_ptr(G, 3, F), F));
        t = gr_is_zero(y, F);
        if (t != T_FALSE)
        {
            flint_printf("FAIL: nonzero (%d), case %d\n", (int) t, which);
            flint_abort();
        }

        gr_heap_clear(x, F);
        gr_heap_clear(y, F);
        gr_vec_clear(G, F);
        gr_ctx_clear(F);
        gr_tower_clear(T);
    }

    gr_ctx_clear(QQ);
    TEST_FUNCTION_END(state);
}

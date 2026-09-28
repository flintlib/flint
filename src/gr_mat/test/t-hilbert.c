/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "arb.h"
#include "gr.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_hilbert, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_ptr t;
        slong R, C, i, j, prec = 0;
        ulong p = 0;
        int which, status, expected;

        R = n_randint(state, 12);
        C = n_randint(state, 12);
        which = n_randint(state, 4);

        if (which == 0)
        {
            gr_ctx_init_fmpq(ctx);
        }
        else if (which == 1)
        {
            /* prime > R + C - 1: all entries exist */
            p = n_nextprime(R + C + n_randint(state, 100), 1);
            GR_MUST_SUCCEED(gr_ctx_init_nmod(ctx, p));
        }
        else if (which == 2)
        {
            /* small prime: 1/(i+j+1) may not exist */
            p = n_nextprime(n_randint(state, 20), 1);
            GR_MUST_SUCCEED(gr_ctx_init_nmod(ctx, p));
        }
        else
        {
            prec = 2 + n_randint(state, 200);
            gr_ctx_init_real_arb(ctx, prec);
        }

        gr_mat_init(A, R, C, ctx);
        t = gr_heap_init(ctx);
        GR_MUST_SUCCEED(gr_mat_randtest(A, state, ctx));

        status = gr_mat_hilbert(A, ctx);

        if (which == 1 || which == 2)
            expected = (R != 0 && C != 0 && (ulong) (R + C - 1) >= p) ? GR_DOMAIN : GR_SUCCESS;
        else
            expected = GR_SUCCESS;

        if (status != expected)
        {
            TEST_FUNCTION_FAIL("status: R = %wd, C = %wd, p = %wu, status = %d, expected = %d\n%{gr_ctx}\n",
                R, C, p, status, expected, ctx);
        }

        if (status == GR_SUCCESS)
        {
            for (i = 0; i < R; i++)
            {
                for (j = 0; j < C; j++)
                {
                    gr_srcptr e = gr_mat_entry_srcptr(A, i, j, ctx);

                    if (which == 3)
                    {
                        /* must agree bitwise with the direct computation,
                           and contain the exact value */
                        fmpq_t q;
                        arb_one(t);
                        arb_div_ui(t, t, i + j + 1, prec);

                        fmpq_init(q);
                        fmpq_set_si(q, 1, i + j + 1);

                        if (!arb_equal(e, t) || !arb_contains_fmpq(e, q))
                        {
                            TEST_FUNCTION_FAIL("arb: R = %wd, C = %wd, i = %wd, j = %wd\n%{gr}\n%{gr}\n",
                                R, C, i, j, e, ctx, t, ctx);
                        }

                        fmpq_clear(q);
                    }
                    else
                    {
                        GR_MUST_SUCCEED(gr_mul_ui(t, e, i + j + 1, ctx));

                        if (gr_is_one(t, ctx) != T_TRUE)
                        {
                            TEST_FUNCTION_FAIL("R = %wd, C = %wd, i = %wd, j = %wd\n%{gr_ctx}\n%{gr}\n",
                                R, C, i, j, ctx, e, ctx);
                        }
                    }
                }
            }
        }

        gr_heap_clear(t, ctx);
        gr_mat_clear(A, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

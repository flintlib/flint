/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"
#include "nfloat.h"

TEST_FUNCTION_START(addmul_submul, state)
{
    gr_ctx_t ctx;
    slong i, len, iter, reps, prec;
    gr_ptr x, y, r1, r2;

    /* test that addmul and submul specializations match the trivial
       algorithm */
    for (prec = FLINT_BITS; prec <= NFLOAT_MAX_LIMBS * FLINT_BITS; prec += FLINT_BITS)
    {
        GR_MUST_SUCCEED(nfloat_ctx_init(ctx, prec, 0));

        reps = (prec <= 256 ? 1000 : 10) * flint_test_multiplier();

        for (iter = 0; iter < reps; iter++)
        {
            len = n_randint(state, 4);

            x = gr_heap_init_vec(len, ctx);
            y = gr_heap_init(ctx);
            r1 = gr_heap_init_vec(len, ctx);
            r2 = gr_heap_init_vec(len, ctx);

            GR_IGNORE(_gr_vec_randtest(x, state, len, ctx));
            GR_IGNORE(gr_randtest(y, state, ctx));
            GR_IGNORE(_gr_vec_randtest(r1, state, len, ctx));
            GR_IGNORE(_gr_vec_set(r2, r1, len, ctx));

            if (n_randint(state, 2))
            {
                GR_MUST_SUCCEED(_gr_vec_addmul_scalar(r1, x, len, y, ctx));

                if (n_randint(state, 2))
                    for (i = 0; i < len; i++)
                        GR_MUST_SUCCEED(gr_addmul(GR_ENTRY(r2, i, ctx->sizeof_elem), 
                                                  GR_ENTRY(x, i, ctx->sizeof_elem), y, ctx));
                else
                    GR_MUST_SUCCEED(gr_generic_vec_scalar_addmul(r2, x, len, y, ctx));
            }
            else
            {
                GR_MUST_SUCCEED(_gr_vec_submul_scalar(r1, x, len, y, ctx));

                if (n_randint(state, 2))
                    for (i = 0; i < len; i++)
                        GR_MUST_SUCCEED(gr_submul(GR_ENTRY(r2, i, ctx->sizeof_elem), 
                                              GR_ENTRY(x, i, ctx->sizeof_elem), y, ctx));
                else
                    GR_MUST_SUCCEED(gr_generic_vec_scalar_submul(r2, x, len, y, ctx));
            }

            if (_gr_vec_equal(r1, r2, len, ctx) != T_TRUE)
            {
                flint_printf("FAIL: %wd\n", prec / FLINT_BITS);
                flint_printf("x = "); _gr_vec_print(x, len, ctx); flint_printf("\n\n");
                flint_printf("y = "); gr_println(y, ctx);
                flint_printf("r1 = "); _gr_vec_print(r1, len, ctx); flint_printf("\n\n");
                flint_printf("r2 = "); _gr_vec_print(r2, len, ctx); flint_printf("\n\n");
                flint_abort();
            }

            gr_heap_clear_vec(x, len, ctx);
            gr_heap_clear(y, ctx);
            gr_heap_clear_vec(r1, len, ctx);
            gr_heap_clear_vec(r2, len, ctx);
        }
        
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

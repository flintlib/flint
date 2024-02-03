/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_get_set_fq_nmod, state)
{
    slong ix, jx;
    int result;

    for (ix = 0; ix < 200 * flint_test_multiplier(); ix++)
    {
        fq_zech_ctx_t ctx;
        fq_zech_t a, b;
        fq_nmod_t c;

        fq_zech_ctx_init_randtest(ctx, state, 3);
        fq_zech_init(a, ctx);
        fq_zech_init(b, ctx);
        fq_nmod_init(c, ctx->fq_nmod_ctx);

        for (jx = 0; jx < 20; jx++)
        {
            fq_zech_randtest(a, state, ctx);
            fq_zech_get_fq_nmod(c, a, ctx);
            fq_zech_set_fq_nmod(b, c, ctx);

            result = (fq_zech_equal(a, b, ctx));
            if (!result)
            {
                flint_printf("FAIL:n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_nmod_print_pretty(c, ctx->fq_nmod_ctx), flint_printf("\n");
                flint_printf("table = %wd\n", ctx->eval_table[a->value]);
                flint_abort();
            }
        }

        fq_zech_clear(a, ctx);
        fq_zech_clear(b, ctx);
        fq_nmod_clear(c, ctx->fq_nmod_ctx);
        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

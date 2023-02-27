/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "ulong_extras.h"
#include "fq_nmod_mpoly.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_is_fq_nmod....");
    fflush(stdout);

    /* Set to random integer and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        fq_nmod_t c, d;
        ulong b;
        slong len;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);

        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(d, ctx->fqctx);

        len = n_randint(state, 10);
        exp_bits = n_randint(state, 200);

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        if (fq_nmod_mpoly_is_fq_nmod(f, ctx))
        {
            fq_nmod_mpoly_get_fq_nmod(c, f, ctx);
            if (!fq_nmod_mpoly_equal_fq_nmod(f, c, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check is_fq_nmod and get_fq_nmod catch constants\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_randtest(c, state, ctx->fqctx);
        fq_nmod_mpoly_set_fq_nmod(f, c, ctx);
        if (!fq_nmod_mpoly_is_fq_nmod(f, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check set_fq_nmod makes is_fq_nmod true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
        fq_nmod_mpoly_get_fq_nmod(d, f, ctx);
        if (!fq_nmod_equal(c, d, ctx->fqctx))
        {
            printf("FAIL\n");
            flint_printf("Check get_fq_nmod matches set_fq_nmod\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        b = n_randlimb(state);
        fq_nmod_set_ui(c, b, ctx->fqctx);
        fq_nmod_mpoly_set_ui(f, b, ctx);
        if (!fq_nmod_mpoly_is_fq_nmod(f, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check set_ui makes is_fq_nmod true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
        fq_nmod_mpoly_get_fq_nmod(d, f, ctx);
        if (!fq_nmod_equal(c, d, ctx->fqctx))
        {
            printf("FAIL\n");
            flint_printf("Check get_fq_nmod matches set_ui\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_clear(d, ctx->fqctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


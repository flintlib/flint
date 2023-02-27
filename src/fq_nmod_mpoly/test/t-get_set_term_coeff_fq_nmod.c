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
#include "fq_nmod_mpoly.h"

int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_term_coeff_fq_nmod....");
    fflush(stdout);

    /* Set coeff and get coeff and compare */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        fq_nmod_t c, d;
        slong len, index;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(d, ctx->fqctx);

        len = n_randint(state, 100) + 1;
        exp_bits = n_randint(state, 200) + 1;
        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
        if (fq_nmod_mpoly_is_zero(f, ctx))
            fq_nmod_mpoly_one(f, ctx);

        for (j = 0; j < 10; j++)
        {
            fq_nmod_randtest(c, state, ctx->fqctx);

            index = n_randint(state, f->length);

            fq_nmod_mpoly_set_term_coeff_fq_nmod(f, index, c, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(d, f, index, ctx);
            if (!fq_nmod_equal(c, d, ctx->fqctx))
            {
                printf("FAIL\n");
                flint_printf("check get and set match\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_clear(d, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

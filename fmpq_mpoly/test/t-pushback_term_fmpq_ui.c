/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("pushback_fmpq_ui....");
    fflush(stdout);

    /* Check pushback matches add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f1, f2, m;
        mp_bitcnt_t coeff_bits, exp_bits;
        ulong * exp;
        slong len;
        fmpq_t c;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        fmpq_mpoly_init(f1, ctx);
        fmpq_mpoly_init(f2, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_init(c);
        exp = (ulong *) flint_malloc(ctx->zctx->minfo->nvars*sizeof(ulong));

        len = n_randint(state, 20);
        coeff_bits = n_randint(state, 200) + 1;
        exp_bits = n_randint(state, FLINT_BITS + 1);

        fmpq_mpoly_zero(f1, ctx);
        fmpq_mpoly_zero(f2, ctx);

        for (j = 0; j < len; j++)
        {
            /* get random term */
            fmpq_randtest(c, state, coeff_bits);
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
                exp[k] = n_randint(state, exp_bits);

            /* add it to f1 */
            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_term_fmpq_ui(m, c, exp, ctx);
            fmpq_mpoly_add(f1, f1, m, ctx);

            /* push it back on f2 */
            fmpq_mpoly_pushback_term_fmpq_ui(f2, c, exp, ctx);
        }

        fmpq_mpoly_sort(f2, ctx);
        fmpq_mpoly_combine_like_terms(f2, ctx);
        fmpq_mpoly_assert_canonical(f2, ctx);

        if (!fmpq_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushback matches add\ni=%wd\n",i,j);
            flint_abort();
        }

        fmpq_clear(c);
        fmpq_mpoly_clear(f1, ctx);
        fmpq_mpoly_clear(f2, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_ctx_clear(ctx);

        flint_free(exp);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


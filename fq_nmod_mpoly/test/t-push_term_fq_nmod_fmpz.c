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
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("push_term_fq_nmod_fmpz....");
    fflush(stdout);

    /* Check pushback matches add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f1, f2, m;
        flint_bitcnt_t exp_bits;
        fmpz ** exp, ** exp2;
        slong len, nvars;
        fq_nmod_t c, c2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f1, ctx);
        fq_nmod_mpoly_init(f2, ctx);
        fq_nmod_mpoly_init(m, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(c2, ctx->fqctx);

        nvars = fq_nmod_mpoly_ctx_nvars(ctx);

        exp = (fmpz **) flint_malloc(nvars*sizeof(fmpz *));
        exp2 = (fmpz **) flint_malloc(nvars*sizeof(fmpz *));
        for (k = 0; k < nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(exp[k]);
            exp2[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(exp2[k]);
        }

        len = n_randint(state, 20);
        exp_bits = n_randint(state, 200);

        fq_nmod_mpoly_zero(f1, ctx);
        fq_nmod_mpoly_zero(f2, ctx);

        for (j = 0; j < len; j++)
        {
            /* get random term */
            fq_nmod_randtest(c, state, ctx->fqctx);
            for (k = 0; k < nvars; k++)
                fmpz_randtest_unsigned(exp[k], state, exp_bits);

            /* add it to f1 */
            fq_nmod_mpoly_zero(m, ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(m, c, exp, ctx);
            fq_nmod_mpoly_add(f1, f1, m, ctx);
            fq_nmod_mpoly_assert_canonical(f1, ctx);

            /* push it back on f2 */
            fq_nmod_mpoly_push_term_fq_nmod_fmpz(f2, c, exp, ctx);

            /* make sure last term matches */
            fq_nmod_mpoly_get_term_coeff_fq_nmod(c2, f2,
                                       fq_nmod_mpoly_length(f2, ctx) - 1, ctx);
            fq_nmod_mpoly_get_term_exp_fmpz(exp2, f2,
                                       fq_nmod_mpoly_length(f2, ctx) - 1, ctx);
            if (!fq_nmod_equal(c, c2, ctx->fqctx))
            {
                printf("FAIL\n");
                flint_printf("Check pushed coefficient matches\ni=%wd, j=%wd\n", i, j);
                flint_abort();
            }
            for (k = 0; k < nvars; k++)
            {
                if (!fmpz_equal(exp[k], exp2[k]))
                {
                    printf("FAIL\n");
                    flint_printf("Check pushed exponent matches\ni=%wd, j=%wd\n", i, j);
                    flint_abort();                    
                }
            }
        }

        fq_nmod_mpoly_sort_terms(f2, ctx);
        fq_nmod_mpoly_combine_like_terms(f2, ctx);
        fq_nmod_mpoly_assert_canonical(f2, ctx);

        if (!fq_nmod_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushed polynomial matches add\ni=%wd\n",i,j);
            flint_abort();
        }

        fq_nmod_clear(c2, ctx->fqctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_clear(f1, ctx);
        fq_nmod_mpoly_clear(f2, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);

        for (k = 0; k < nvars; k++)
        {
            fmpz_clear(exp2[k]);
            flint_free(exp2[k]); 
            fmpz_clear(exp[k]);
            flint_free(exp[k]); 
        }
        flint_free(exp2);
        flint_free(exp);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

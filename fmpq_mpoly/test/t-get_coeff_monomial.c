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
#include "ulong_extras.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("get_coeff_monomial....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, m;
        fmpz ** exp;
        fmpq_t cm, ce, q;
        slong len;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_init(cm);
        fmpq_init(ce);
        fmpq_init(q);
        fmpq_one(q); /* anything nonzero is ok */

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);
        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        /* check a random monomial - this also randomizes m->bits */
        exp = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz *));
        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(exp[k]);
            fmpz_randtest_unsigned(exp[k], state, exp_bits1);
        }
        fmpq_mpoly_zero(m, ctx);
        fmpq_mpoly_set_term_fmpq_fmpz(m, q, exp, ctx);
        fmpq_mpoly_get_coeff_monomial(cm, f, m, ctx);
        fmpq_mpoly_get_term_fmpq_fmpz(ce, f, exp, ctx);
        if (!fmpq_equal(cm, ce))
        {
            flint_printf("FAIL\ni = %wd\n", i);
            flint_abort();
        }

        /* check all monomials in f */
        for (j = 0; j < fmpq_mpoly_length(f, ctx); j++)
        {
            fmpq_mpoly_get_monomial_fmpz(exp, f, j, ctx);

            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_term_fmpq_fmpz(m, q, exp, ctx);
            fmpq_mpoly_get_coeff_monomial(cm, f, m, ctx);
            fmpq_mpoly_get_term_fmpq_fmpz(ce, f, exp, ctx);
            if (!fmpq_equal(cm, ce))
            {
                flint_printf("FAIL\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_clear(exp[k]);
            flint_free(exp[k]); 
        }
        flint_free(exp);

        fmpq_clear(q);
        fmpq_clear(cm);
        fmpq_clear(ce);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


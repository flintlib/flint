/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j, k;
    int result;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_term_exp_fmpz....");
    fflush(stdout);

    /* check get and set match */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        slong len, nvars, index;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        nvars = fmpq_mpoly_ctx_nvars(ctx);

        do {
            len = n_randint(state, 50);
            exp_bits = n_randint(state, 100) + 1;
            coeff_bits = n_randint(state, 100);

            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        } while (fmpq_mpoly_length(f, ctx) == 0);

        for (j = 0; j < 10; j++)
        {
            fmpz ** exp1 = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz*));
            fmpz ** exp2 = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz*));

            for (k = 0; k < nvars; k++)
            {
                exp1[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
                exp2[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
                fmpz_init(exp1[k]);
                fmpz_init(exp2[k]);
                fmpz_randtest_unsigned(exp1[k], state, 200);
            }

            index = n_randint(state, fmpq_mpoly_length(f, ctx));

            fmpq_mpoly_set_term_exp_fmpz(f, index, exp1, ctx);

            if (!mpoly_monomials_valid_test(f->zpoly->exps, f->zpoly->length, f->zpoly->bits, ctx->zctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents invalid");

            if (mpoly_monomials_overflow_test(f->zpoly->exps, f->zpoly->length, f->zpoly->bits, ctx->zctx->minfo))
                flint_throw(FLINT_ERROR, "Polynomial exponents overflow");

            fmpq_mpoly_get_term_exp_fmpz(exp2, f, index, ctx);

            result = 1;
            for (k = 0; k < nvars; k++)
            {
                result = result && fmpz_equal(exp1[k], exp2[k]);
                fmpz_clear(exp1[k]);
                fmpz_clear(exp2[k]);
                flint_free(exp1[k]); 
                flint_free(exp2[k]); 
            }

            if (!result)
            {
                flint_printf("FAIL\ncheck get and set match\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }

            flint_free(exp1);
            flint_free(exp2);
       }

       fmpq_mpoly_clear(f, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

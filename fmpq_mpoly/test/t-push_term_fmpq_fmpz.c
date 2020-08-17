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

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("push_back_fmpq_fmpz....");
    fflush(stdout);

    /* Check pushback matches add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f1, f2, m;
        flint_bitcnt_t coeff_bits, exp_bits;
        fmpz ** exp, ** exp2;
        slong len, nvars;
        fmpq_t c, c2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        fmpq_mpoly_init(f1, ctx);
        fmpq_mpoly_init(f2, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_init(c);
        fmpq_init(c2);

        nvars = fmpq_mpoly_ctx_nvars(ctx);

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
        coeff_bits = n_randint(state, 100) + 1;
        exp_bits = n_randint(state, 200);

        fmpq_mpoly_zero(f1, ctx);
        fmpq_mpoly_zero(f2, ctx);

        for (j = 0; j < len; j++)
        {
            /* get random term */
            fmpq_randtest(c, state, coeff_bits + 1);
            for (k = 0; k < nvars; k++)
                fmpz_randtest_unsigned(exp[k], state, exp_bits);

            /* add it to f1 */
            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_coeff_fmpq_fmpz(m, c, exp, ctx);
            fmpq_mpoly_add(f1, f1, m, ctx);
            fmpq_mpoly_assert_canonical(f1, ctx);

            /* push it back on f2 */
            if (fmpz_is_one(fmpq_denref(c)))
                fmpq_mpoly_push_term_fmpz_fmpz(f2, fmpq_numref(c), exp, ctx);
            else
                fmpq_mpoly_push_term_fmpq_fmpz(f2, c, exp, ctx);

            /* make sure last term matches */
            fmpq_mpoly_get_term_coeff_fmpq(c2, f2, fmpq_mpoly_length(f2, ctx) - 1, ctx);
            fmpq_mpoly_get_term_exp_fmpz(exp2, f2, fmpq_mpoly_length(f2, ctx) - 1, ctx);
            if (!fmpq_equal(c, c2))
            {
                printf("FAIL\n");
                flint_printf("Check pushed coefficient matches\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
            for (k = 0; k < nvars; k++)
            {
                if (!fmpz_equal(exp[k], exp2[k]))
                {
                    printf("FAIL\n");
                    flint_printf("Check pushed exponent matches\ni = %wd, j = %wd\n", i, j);
                    flint_abort();                    
                }
            }
        }

        fmpq_mpoly_sort_terms(f2, ctx);
        fmpq_mpoly_combine_like_terms(f2, ctx);
        fmpq_mpoly_assert_canonical(f2, ctx);

        if (!fmpq_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushed polynomial matches add\ni = %wd\n",i);
            flint_abort();
        }

        fmpq_clear(c2);
        fmpq_clear(c);
        fmpq_mpoly_clear(f1, ctx);
        fmpq_mpoly_clear(f2, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_ctx_clear(ctx);

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


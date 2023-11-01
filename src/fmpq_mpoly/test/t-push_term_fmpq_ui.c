/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_push_term_fmpq_ui, state)
{
    slong i, j, k;

    /* Check pushback matches add */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f1, f2, m;
        flint_bitcnt_t coeff_bits, exp_bits;
        ulong * exp, * exp2;
        slong len, nvars;
        fmpq_t c, c2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        fmpq_mpoly_init(f1, ctx);
        fmpq_mpoly_init(f2, ctx);
        fmpq_mpoly_init(m, ctx);
        fmpq_init(c);
        fmpq_init(c2);

        nvars = fmpq_mpoly_ctx_nvars(ctx);

        exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
        exp2 = (ulong *) flint_malloc(nvars*sizeof(ulong));

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
                exp[k] = n_randint(state, exp_bits);

            /* add it to f1 */
            fmpq_mpoly_zero(m, ctx);
            fmpq_mpoly_set_coeff_fmpq_ui(m, c, exp, ctx);
            fmpq_mpoly_add(f1, f1, m, ctx);
            fmpq_mpoly_assert_canonical(f1, ctx);

            /* push it back on f2 */
            if (fmpz_is_one(fmpq_denref(c)))
                fmpq_mpoly_push_term_fmpz_ui(f2, fmpq_numref(c), exp, ctx);
            else
                fmpq_mpoly_push_term_fmpq_ui(f2, c, exp, ctx);

            /* make sure last term matches */
            fmpq_mpoly_get_term_coeff_fmpq(c2, f2, fmpq_mpoly_length(f2, ctx) - 1, ctx);
            fmpq_mpoly_get_term_exp_ui(exp2, f2, fmpq_mpoly_length(f2, ctx) - 1, ctx);
            if (!fmpq_equal(c, c2))
            {
                printf("FAIL\n");
                flint_printf("Check pushed coefficient matches\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
            for (k = 0; k < nvars; k++)
            {
                if (exp[k] != exp2[k])
                {
                    printf("FAIL\n");
                    flint_printf("Check pushed exponent matches\ni = %wd, j = %wd\n", i, j);
                    fflush(stdout);
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
            flint_printf("Check pushed polynomial matches add\ni=%wd\n",i,j);
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(c2);
        fmpq_clear(c);
        fmpq_mpoly_clear(f1, ctx);
        fmpq_mpoly_clear(f2, ctx);
        fmpq_mpoly_clear(m, ctx);
        fmpq_mpoly_ctx_clear(ctx);

        flint_free(exp2);
        flint_free(exp);
    }

    TEST_FUNCTION_END(state);
}

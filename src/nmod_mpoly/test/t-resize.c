/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_resize, state)
{
    slong i, j, k;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f1, f2, f3, m;
        flint_bitcnt_t exp_bits;
        ulong * exp, * exp2;
        slong len, nvars;
        mp_limb_t c;
        mp_limb_t modulus;

        modulus = n_randtest_bits(state, n_randint(state, FLINT_BITS) + 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f1, ctx);
        nmod_mpoly_init(f2, ctx);
        nmod_mpoly_init(f3, ctx);
        nmod_mpoly_init(m, ctx);

        nvars = nmod_mpoly_ctx_nvars(ctx);

        exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
        exp2 = (ulong *) flint_malloc(nvars*sizeof(ulong));

        len = n_randint(state, 200);
        exp_bits = n_randint(state, FLINT_BITS) + 1;

        nmod_mpoly_zero(f1, ctx);
        nmod_mpoly_zero(f2, ctx);
        nmod_mpoly_zero(f3, ctx);

        for (j = 0; j < len; j++)
        {
            /* get random term */
            c = n_randlimb(state);
            for (k = 0; k < nvars; k++)
                exp[k] = n_randtest_bits(state, n_randint(state, exp_bits) + 1);

            /* add it to f1 */
            nmod_mpoly_zero(m, ctx);
            nmod_mpoly_set_coeff_ui_ui(m, c, exp, ctx);
            nmod_mpoly_add(f1, f1, m, ctx);

            /* push it back on f2 */
            nmod_mpoly_push_term_ui_ui(f2, c, exp, ctx);

            /* manually push it on f3 */
            nmod_mpoly_resize(f3, j + 1 + n_randint(state, 10), ctx);
            nmod_mpoly_set_term_coeff_ui(f3, j, c, ctx);
            nmod_mpoly_set_term_exp_ui(f3, j, exp, ctx);
        }

        nmod_mpoly_sort_terms(f2, ctx);
        nmod_mpoly_combine_like_terms(f2, ctx);
        nmod_mpoly_assert_canonical(f2, ctx);

        nmod_mpoly_sort_terms(f3, ctx);
        nmod_mpoly_combine_like_terms(f3, ctx);
        nmod_mpoly_assert_canonical(f3, ctx);

        if (!nmod_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushback matches add\ni=%wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        if (!nmod_mpoly_equal(f1, f3, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check resize+setterm matches add\ni=%wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f1, ctx);
        nmod_mpoly_clear(f2, ctx);
        nmod_mpoly_clear(f3, ctx);
        nmod_mpoly_clear(m, ctx);
        nmod_mpoly_ctx_clear(ctx);

        flint_free(exp2);
        flint_free(exp);
    }

    TEST_FUNCTION_END(state);
}

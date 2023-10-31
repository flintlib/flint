/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_get_set_term_coeff_ui, state)
{
    slong i, j;

    /* Set coeff and get coeff and compare */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f;
        slong len, index;
        flint_bitcnt_t exp_bits;
        mp_limb_t c, d;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;
        nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        if (f->length > 0)
        {
            for (j = 0; j < 10; j++)
            {
                c = n_randint(state, modulus);

                index = n_randint(state, f->length);

                nmod_mpoly_set_term_coeff_ui(f, index, c, ctx);
                d = nmod_mpoly_get_term_coeff_ui(f, index, ctx);
                if (c != d)
                {
                    printf("FAIL\n");
                    flint_printf("check get and set match\ni = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                if (*nmod_mpoly_term_coeff_ref(f, index, ctx) != d)
                {
                    printf("FAIL\n");
                    flint_printf("check reference match\ni = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

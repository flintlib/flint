/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, k, result;
    FLINT_TEST_INIT(state);

    flint_printf("get/set_term....");
    fflush(stdout);

    /* Check _ui_fmpz */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f;
        mp_bitcnt_t exp_bits;
        slong len;
        mp_limb_t c, d;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;

        nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz * exp = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));

            c = n_randint(state, modulus);
            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_init(exp + k);
                fmpz_randtest_unsigned(exp + k, state, 200);
            }

            _nmod_mpoly_set_term_ui_fmpz(f, c, exp, ctx);
            nmod_mpoly_test(f, ctx);
            d = _nmod_mpoly_get_term_ui_fmpz(f, exp, ctx);
            result = (c == d);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _ui_fmpz\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            for (k = 0; k < ctx->minfo->nvars; k++)
                fmpz_clear(exp + k);

            flint_free(exp);
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check _ui_ui */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f;
        mp_bitcnt_t exp_bits;
        slong len;
        mp_limb_t c, d;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;

        nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            ulong * exp = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

            c = n_randint(state, modulus);
            for (k = 0; k < ctx->minfo->nvars; k++)
                exp[k] = n_randtest(state);

            nmod_mpoly_set_term_ui_ui(f, c, exp, ctx);
            nmod_mpoly_test(f, ctx);
            d = nmod_mpoly_get_term_ui_ui(f, exp, ctx);
            result = (c == d);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check _ui_ui\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            flint_free(exp);
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


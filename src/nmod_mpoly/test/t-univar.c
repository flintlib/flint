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

TEST_FUNCTION_START(nmod_mpoly_univar, state)
{
    slong i, j, k;

    /* Check mpoly -> mpoly_univar -> mpoly */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        nmod_mpoly_univar_t fx, gx;
        slong len1, len2, n;
        flint_bitcnt_t exp_bits1, exp_bits2, bits;
        mp_limb_t modulus;

        modulus = n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 2;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_univar_init(fx, ctx);
        nmod_mpoly_univar_init(gx, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        exp_bits1 = n_randint(state, 3*FLINT_BITS) + 1;
        exp_bits2 = n_randint(state, 3*FLINT_BITS) + 1;

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            nmod_mpoly_to_univar(fx, f, j, ctx);
            nmod_mpoly_univar_assert_canonical(fx, ctx);
            nmod_mpoly_from_univar(g, fx, j, ctx);
            nmod_mpoly_assert_canonical(g, ctx);

            if (!nmod_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check mpoly -> mpoly_univar -> mpoly\ni: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }

            bits = mpoly_fix_bits(f->bits + n_randint(state, FLINT_BITS), ctx->minfo);
            _nmod_mpoly_from_univar(h, bits, fx, j, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            if (h->bits != bits || !nmod_mpoly_equal(f, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check mpoly -> mpoly_univar -> mpoly with bits\ni: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }

            if (!nmod_mpoly_univar_degree_fits_si(fx, ctx))
                continue;

            n = nmod_mpoly_univar_length(fx, ctx);
            nmod_mpoly_univar_fit_length(gx, n, ctx);
            gx->length = n;
            for (k = 0; k < n; k++)
            {
                nmod_mpoly_univar_swap_term_coeff(gx->coeffs + k, fx, k, ctx);
                fmpz_set_si(gx->exps + k, nmod_mpoly_univar_get_term_exp_si(fx, k, ctx));
            }

            nmod_mpoly_from_univar(g, gx, j, ctx);
            if (!nmod_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check construction\ni: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_univar_clear(gx, ctx);
        nmod_mpoly_univar_clear(fx, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

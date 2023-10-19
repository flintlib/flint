/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_univar, state)
{
    slong i, j, k;

    /* Check mpoly -> mpoly_univar -> mpoly */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        fmpz_mpoly_univar_t fx, gx;
        slong len1, len2, n;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2, bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 15);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_univar_init(fx, ctx);
        fmpz_mpoly_univar_init(gx, ctx);

        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        exp_bits1 = n_randint(state, (i%2) ? 4 : 3*FLINT_BITS) + 1;
        exp_bits2 = n_randint(state, (i%2) ? 4 : 3*FLINT_BITS) + 1;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_to_univar(fx, f, j, ctx);
            fmpz_mpoly_univar_assert_canonical(fx, ctx);
            fmpz_mpoly_from_univar(g, fx, j, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!fmpz_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL: Check mpoly -> mpoly_univar -> mpoly\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }

            bits = mpoly_fix_bits(f->bits + n_randint(state, FLINT_BITS), ctx->minfo);
            _fmpz_mpoly_from_univar(h, bits, fx, j, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            if (h->bits != bits || !fmpz_mpoly_equal(f, h, ctx))
            {
                flint_printf("FAIL: Check mpoly -> mpoly_univar -> mpoly with bits\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }

            if (!fmpz_mpoly_univar_degree_fits_si(fx, ctx))
                continue;

            n = fmpz_mpoly_univar_length(fx, ctx);
            fmpz_mpoly_univar_fit_length(gx, n, ctx);
            gx->length = n;
            for (k = 0; k < n; k++)
            {
                fmpz_mpoly_univar_swap_term_coeff(gx->coeffs + k, fx, k, ctx);
                fmpz_set_si(gx->exps + k, fmpz_mpoly_univar_get_term_exp_si(fx, k, ctx));
            }

            fmpz_mpoly_from_univar(g, gx, j, ctx);
            if (!fmpz_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL: Check construction\n");
                flint_printf("i: %wd  j: %wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_univar_clear(gx, ctx);
        fmpz_mpoly_univar_clear(fx, ctx);

        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

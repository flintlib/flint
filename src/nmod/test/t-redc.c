/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"

TEST_FUNCTION_START(nmod_redc, state)
{
    int i;

    for (i = 0; i < 100000 * flint_test_multiplier(); i++)
    {
        nmod_t mod;
        nmod_redc_ctx_t ctx;
        ulong n, x, y, z, x_red, y_red, z_red, z_back, e;

        n = n_randtest(state) | 1;

        nmod_init(&mod, n);

        if (n_randint(state, 2))
            nmod_redc_ctx_init_nmod(ctx, mod);
        else
            nmod_redc_ctx_init_ui(ctx, n);

        x = n_randtest(state) % n;
        y = n_randtest(state) % n;
        e = n_randtest_not_zero(state);

        x_red = nmod_redc_set_nmod(x, ctx);
        y_red = nmod_redc_set_nmod(y, ctx);

        z = nmod_add(x, y, mod);
        z_red = nmod_redc_add(x_red, y_red, ctx);
        z_back = nmod_redc_get_nmod(z_red, ctx);

        if (z_back != z)
        {
            TEST_FUNCTION_FAIL("nmod_redc_add\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
        }

        z = nmod_mul(x, y, mod);
        z_red = nmod_redc_mul(x_red, y_red, ctx);
        z_back = nmod_redc_get_nmod(z_red, ctx);

        if (z_back != z)
        {
            TEST_FUNCTION_FAIL("nmod_redc_mul\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
        }

        z = _nmod_pow_ui_binexp(x, e, mod);
        z_red = _nmod_redc_pow_ui(x_red, e, ctx);
        z_back = nmod_redc_get_nmod(z_red, ctx);

        if (z_back != z)
        {
            TEST_FUNCTION_FAIL("_nmod_redc_pow_ui\nn = %wu\nnred = %wu\nx = %wu\ne = %wu\nz = %wu\nx_red = %wu\nz_red = %wu\nz_back = %wu\n",
                n, ctx->nred, x, e, z, x_red, z_red, z_back);
        }

        z = (e == 0) ? nmod_set_ui(1, mod) : _nmod_pow_ui_binexp(2 % mod.n, e, mod);
        z_red = _nmod_redc_2_pow_ui(e, ctx);
        z_back = nmod_redc_get_nmod(z_red, ctx);

        if (z_back != z)
        {
            TEST_FUNCTION_FAIL("_nmod_redc_2_pow_ui\nn = %wu\nnred = %wu\ne = %wu\nz = %wu\nz_red = %wu\nz_back = %wu\n",
                n, ctx->nred, e, z, z_red, z_back);
        }

        if (nmod_redc_can_use_fast(ctx))
        {
            if (n_randint(state, 2)) x_red += n;
            if (n_randint(state, 2)) y_red += n;

            z = nmod_add(x, y, mod);
            z_red = nmod_redc_fast_add(x_red, y_red, ctx);
            z_back = nmod_redc_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_fast_add\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
            }

            z = nmod_mul(x, y, mod);
            z_red = nmod_redc_fast_mul(x_red, y_red, ctx);
            z_back = nmod_redc_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_fast_mul\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
            }

            z = nmod_add(x, x, mod);
            z_red = nmod_redc_fast_mul_two(x_red, ctx);
            z_back = nmod_redc_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_fast_mul_two\nn = %wu\nnred = %wu\nx = %wu\nz = %wu\nx_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, z, x_red, z_red, z_back);
            }

            z = _nmod_pow_ui_binexp(x, e, mod);
            z_red = _nmod_redc_fast_pow_ui(x_red, e, ctx);
            z_back = nmod_redc_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("_nmod_redc_fast_pow_ui\nn = %wu\nnred = %wu\nx = %wu\ne = %wu\nz = %wu\nx_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, e, z, x_red, z_red, z_back);
            }

            z = (e == 0) ? nmod_set_ui(1, mod) : _nmod_pow_ui_binexp(2 % mod.n, e, mod);
            z_red = _nmod_redc_fast_2_pow_ui(e, ctx);
            z_back = nmod_redc_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("_nmod_redc_fast_2_pow_ui\nn = %wu\nnred = %wu\ne = %wu\nz = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, e, z, z_red, z_back);
            }
        }

        if (n < UWORD(1) << (FLINT_BITS / 2 - 1))
        {
            if (n_randint(state, 2))
                nmod_redc_half_ctx_init_nmod(ctx, mod);
            else
                nmod_redc_half_ctx_init_ui(ctx, n);

            x = n_randtest(state) % n;
            y = n_randtest(state) % n;
            e = n_randtest_not_zero(state);

            x_red = nmod_redc_half_set_nmod(x, ctx);
            y_red = nmod_redc_half_set_nmod(y, ctx);

            z = nmod_add(x, y, mod);
            z_red = nmod_redc_half_add(x_red, y_red, ctx);
            z_back = nmod_redc_half_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_half_add\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
            }

            z = nmod_sub(x, y, mod);
            z_red = nmod_redc_half_sub(x_red, y_red, ctx);
            z_back = nmod_redc_half_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_half_sub\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
            }

            z = nmod_mul(x, y, mod);
            z_red = nmod_redc_half_mul(x_red, y_red, ctx);
            z_back = nmod_redc_half_get_nmod(z_red, ctx);

            if (z_back != z)
            {
                TEST_FUNCTION_FAIL("nmod_redc_half_mul\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                    n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
            }

            if (n < UWORD(1) << (FLINT_BITS / 2 - 2))
            {
                if (n_randint(state, 2)) x_red += n;
                if (n_randint(state, 2)) y_red += n;

                z = nmod_add(x, y, mod);
                z_red = nmod_redc_half_fast_add(x_red, y_red, ctx);
                z_back = nmod_redc_half_get_nmod(z_red, ctx);

                if (z_back != z)
                {
                    TEST_FUNCTION_FAIL("nmod_redc_half_fast_add\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                        n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
                }

                z = nmod_sub(x, y, mod);
                z_red = nmod_redc_half_fast_sub(x_red, y_red, ctx);
                z_back = nmod_redc_half_get_nmod(z_red, ctx);

                if (z_back != z)
                {
                    TEST_FUNCTION_FAIL("nmod_redc_half_fast_sub\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                        n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
                }

                z = nmod_mul(x, y, mod);
                z_red = nmod_redc_half_fast_mul(x_red, y_red, ctx);
                z_back = nmod_redc_half_get_nmod(z_red, ctx);

                if (z_back != z)
                {
                    TEST_FUNCTION_FAIL("nmod_redc_half_fast_mul\nn = %wu\nnred = %wu\nx = %wu\ny = %wu\nz = %wu\nx_red = %wu\ny_red = %wu\nz_red = %wu\nz_back = %wu\n",
                        n, ctx->nred, x, y, z, x_red, y_red, z_red, z_back);
                }
            }
        }
    }

    TEST_FUNCTION_END(state);
}

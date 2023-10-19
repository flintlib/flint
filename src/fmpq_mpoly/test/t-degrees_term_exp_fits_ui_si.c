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
#include "ulong_extras.h"

TEST_FUNCTION_START(fmpq_mpoly_degrees_term_exp_fits_ui_si, state)
{
    slong i, j, k;
    int result;

    /* basic corner cases */
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        const char * vars[] = {"x","y","z"};

        fmpq_mpoly_ctx_init(ctx, 3, ORD_DEGLEX);
        fmpq_mpoly_init(f, ctx);

        if (FLINT_BITS == 64)
        {
            fmpq_mpoly_set_str_pretty(f, "x^9223372036854775807", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 1\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^9223372036854775807", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 2\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^9223372036854775808", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 3\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^9223372036854775808", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 4\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^18446744073709551615", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^18446744073709551615", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^18446744073709551616", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^18446744073709551616", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                fflush(stdout);
                flint_abort();
            }

        } else if (FLINT_BITS == 32)
        {
            fmpq_mpoly_set_str_pretty(f, "x^2147483647", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 1\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483647*y^2147483647", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 2\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483648", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 3\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483648*y^2147483648", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 4\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967295", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967295*y^4294967295", vars, ctx);
            if (!fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967296", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967296*y^4294967296", vars, ctx);
            if (fmpq_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                fflush(stdout);
                flint_abort();
            }

        } else
        {
            printf("FAIL\nFLINT_BITS is not 64 or 32\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check fit */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpz ** exp, ** deg;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        exp = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz *));
        deg = (fmpz **) flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz *));
        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            deg[k] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(exp[k]);
            fmpz_init(deg[k]);
            fmpz_set_si(deg[k], -WORD(1));
        }

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        coeff_bits = n_randint(state, 30);

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < fmpq_mpoly_length(f, ctx); j++)
        {
            fmpq_mpoly_get_term_exp_fmpz(exp, f, j, ctx);
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
                if (fmpz_cmp(deg[k], exp[k]) < 0)
                    fmpz_set(deg[k], exp[k]);

            result = 1;
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
                result = result && fmpz_fits_si(exp[k]);
            if (result != fmpq_mpoly_term_exp_fits_si(f, j, ctx))
            {
                flint_printf("FAIL\nCheck monomial_fit_si\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            result = 1;
            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
                result = result && fmpz_abs_fits_ui(exp[k]);
            if (result != fmpq_mpoly_term_exp_fits_ui(f, j, ctx))
            {
                flint_printf("FAIL\nCheck monomial_fit_ui\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        result = 1;
        fmpq_mpoly_degrees_fmpz(exp, f, ctx);
        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            if (!fmpz_equal(exp[k], deg[k]))
            {
                flint_printf("FAIL\nCheck degree computation\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
            result = result && fmpz_fits_si(exp[k]);
        }
        if (result != fmpq_mpoly_degrees_fit_si(f, ctx))
        {
            flint_printf("FAIL\nCheck degrees_fit_si\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_clear(deg[k]);
            flint_free(deg[k]);
            fmpz_clear(exp[k]);
            flint_free(exp[k]);
        }
        flint_free(deg);
        flint_free(exp);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

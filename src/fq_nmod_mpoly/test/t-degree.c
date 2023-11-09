/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

/* check various degree operations */
void _check_degrees(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    int fits_si;
    fmpz ** degs;
    slong * degs_si;
    slong nvars = ctx->minfo->nvars;

    degs = (fmpz **) flint_malloc(nvars*sizeof(fmpz *));
    for (i = 0; i < nvars; i++)
    {
        degs[i] = flint_malloc(sizeof(fmpz));
        fmpz_init(degs[i]);
    }

    fq_nmod_mpoly_degrees_fmpz(degs, A, ctx);

    fits_si = 1;
    for (i = 0; i < nvars; i++)
    {
        fits_si = fits_si && fmpz_fits_si(degs[i]);
    }

    if (fits_si != fq_nmod_mpoly_degrees_fit_si(A, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check degrees_fit_si\n");
        fflush(stdout);
        flint_abort();
    }

    if (fits_si)
    {
        degs_si = (slong *) flint_malloc(nvars*sizeof(slong));
        fq_nmod_mpoly_degrees_si(degs_si, A, ctx);
        for (i = 0; i < nvars; i++)
        {
            if (degs_si[i] != fmpz_get_si(degs[i]))
            {
                printf("FAIL\n");
                flint_printf("Check degrees_si\n");
                fflush(stdout);
                flint_abort();
            }
            if (degs_si[i] != fq_nmod_mpoly_degree_si(A, i, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check individual degree_si\n");
                fflush(stdout);
                flint_abort();
            }
        }
        flint_free(degs_si);
    }

    for (i = 0; i < nvars; i++)
    {
        fmpz_t degi;
        fmpz_init(degi);

        fq_nmod_mpoly_degree_fmpz(degi, A, i, ctx);

        if (!fmpz_equal(degi, degs[i]))
        {
            printf("FAIL\n");
            flint_printf("Check individual degree_fmpz\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_clear(degi);
    }

    for (i = 0; i < nvars; i++)
    {
        fmpz_clear(degs[i]);
        flint_free(degs[i]);
    }
    flint_free(degs);
}

TEST_FUNCTION_START(fq_nmod_mpoly_degree, state)
{
    int i, j;

    /* Check degree does not go up under addition */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_add(h, f, g, ctx);

            _check_degrees(h, ctx);
            _check_degrees(f, ctx);
            _check_degrees(g, ctx);

            fq_nmod_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fq_nmod_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fq_nmod_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if ((fmpz_cmp(hdeg, fdeg) > 0) && (fmpz_cmp(hdeg, gdeg) > 0))
            {
                printf("FAIL\n");
                flint_printf("Check degree does not go up under addition\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check degree adds under multiplication */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int ok;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_mul(h, f, g, ctx);

            _check_degrees(h, ctx);
            _check_degrees(f, ctx);
            _check_degrees(g, ctx);

            fq_nmod_mpoly_degree_fmpz(hdeg, h, j, ctx);
            fq_nmod_mpoly_degree_fmpz(fdeg, f, j, ctx);
            fq_nmod_mpoly_degree_fmpz(gdeg, g, j, ctx);

            if (fq_nmod_mpoly_is_zero(g, ctx) || fq_nmod_mpoly_is_zero(f, ctx))
            {
                ok = (fmpz_cmp_si(hdeg, -WORD(1)) == 0);
            }
            else
            {
                fmpz_sub(hdeg, hdeg, fdeg);
                ok = (fmpz_cmp(hdeg, gdeg) == 0);
            }

            if (!ok)
            {
                printf("FAIL\n");
                flint_printf("Check degree adds under multiplication\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

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
#include "fmpq_mpoly.h"


/* foolproof way to check totdeg_check is correct */
void _check_total_degree(const fmpz_t totdeg_check, const fmpq_mpoly_t A,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    slong i, j, N;
    fmpz_t tot, totdeg;
    fmpz * tmp_exps;
    const ulong * exps = A->zpoly->exps;
    slong len = A->zpoly->length;
    mp_bitcnt_t bits = A->zpoly->bits;
    const mpoly_ctx_struct * mctx = ctx->zctx->minfo;
    TMP_INIT;

    fmpz_init(totdeg);
    fmpz_set_si(totdeg, -WORD(1));

    TMP_START;
    fmpz_init(tot);
    tmp_exps = (fmpz *) TMP_ALLOC(mctx->nvars*sizeof(fmpz));
    for (j = 0; j < mctx->nvars; j++)
        fmpz_init(tmp_exps + j);
    N = mpoly_words_per_exp(bits, mctx);
    for (i = 0; i < len; i++)
    {
        mpoly_get_monomial_ffmpz(tmp_exps, exps + N*i, bits, mctx);
        fmpz_zero(tot);
        for (j = 0; j < mctx->nvars; j++)
            fmpz_add(tot, tot, tmp_exps + j);
        if (fmpz_cmp(totdeg, tot) < 0)
            fmpz_swap(totdeg, tot);
    }
    fmpz_clear(tot);
    for (j = 0; j < mctx->nvars; j++)
        fmpz_clear(tmp_exps + j);

    if (!fmpz_equal(totdeg_check, totdeg))
        flint_throw(FLINT_ERROR, "Total degree is wrong");

    fmpz_clear(totdeg);

    TMP_END;
}


int
main(void)
{
    int i, j;

    FLINT_TEST_INIT(state);

    flint_printf("total_degree....");
    fflush(stdout);

    /* Check total_degree does not go up under addition */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_add(h, f, g, ctx);

            fmpq_mpoly_total_degree_fmpz(hdeg, h, ctx);
            fmpq_mpoly_total_degree_fmpz(fdeg, f, ctx);
            fmpq_mpoly_total_degree_fmpz(gdeg, g, ctx);
            _check_total_degree(hdeg, h, ctx);
            _check_total_degree(fdeg, f, ctx);
            _check_total_degree(gdeg, g, ctx);

            if ((fmpz_cmp(hdeg, fdeg) > 0) && (fmpz_cmp(hdeg, gdeg) > 0))
            {
                printf("FAIL\n");
                flint_printf("Check degree does not go up under addition\ni: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check total_degree adds under multiplication */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpz_t fdeg, gdeg, hdeg;
        slong len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpz_init(fdeg);
        fmpz_init(gdeg);
        fmpz_init(hdeg);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_mul(h, f, g, ctx);

            fmpq_mpoly_total_degree_fmpz(hdeg, h, ctx);
            fmpq_mpoly_total_degree_fmpz(fdeg, f, ctx);
            fmpq_mpoly_total_degree_fmpz(gdeg, g, ctx);
            _check_total_degree(hdeg, h, ctx);
            _check_total_degree(fdeg, f, ctx);
            _check_total_degree(gdeg, g, ctx);

            if (fmpq_mpoly_is_zero(f, ctx) || fmpq_mpoly_is_zero(g, ctx))
            {
                ok = fmpz_equal_si(hdeg, -WORD(1))
                     && (fmpz_equal_si(fdeg, -WORD(1))
                         || fmpz_equal_si(gdeg, -WORD(1)));
            }
            else
            {
                fmpz_add(gdeg, gdeg, fdeg);
                ok = fmpz_equal(hdeg, gdeg);
            }

            if (!ok)
            {
                printf("FAIL\n");
                flint_printf("Check degree adds under multiplication\ni: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_clear(fdeg);
        fmpz_clear(gdeg);
        fmpz_clear(hdeg);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}


/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_mul_fmpz, state)
{
    int j, i, result;
    fq_zech_ctx_t ctx;

    for (j = 0; j < 50; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        for (i = 0; i < 200; i++)
        {
            fmpz_t x;
            fq_nmod_t aa, bb;
            fq_zech_t a, b, c;

            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_nmod_init(bb, ctx->fq_nmod_ctx);

            fmpz_init(x);
            fmpz_randtest_mod_signed(x, state, fq_zech_ctx_prime(ctx));

            fq_nmod_randtest(aa, state, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(a, aa, ctx);

            fq_nmod_mul_fmpz(bb, aa, x, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(b, bb, ctx);

            fq_zech_mul_fmpz(c, a, x, ctx);

            result = (fq_zech_equal(b, c, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("aa = ");
                fq_nmod_print_pretty(aa, ctx->fq_nmod_ctx);
                flint_printf("\n");
                flint_printf("a = ");
                fq_zech_print_pretty(a, ctx);
                flint_printf("\n");
                flint_printf("b = ");
                fq_zech_print_pretty(b, ctx);
                flint_printf("\n");
                flint_printf("c = ");
                fq_zech_print_pretty(c, ctx);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(x);
            fq_nmod_clear(bb, ctx->fq_nmod_ctx);
            fq_nmod_clear(aa, ctx->fq_nmod_ctx);
        }

        for (i = 0; i < 200; i++)
        {
            fmpz_t x;
            fq_nmod_t aa, bb;
            fq_zech_t a, b;

            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_nmod_init(bb, ctx->fq_nmod_ctx);

            fmpz_init(x);
            fmpz_randtest_mod_signed(x, state, fq_zech_ctx_prime(ctx));

            fq_nmod_randtest(aa, state, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(a, aa, ctx);

            fq_nmod_mul_fmpz(bb, aa, x, ctx->fq_nmod_ctx);
            fq_zech_set_fq_nmod(b, bb, ctx);

            fq_zech_mul_fmpz(a, a, x, ctx);

            result = (fq_zech_equal(b, a, ctx));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                fq_zech_ctx_print(ctx);
                flint_printf("\n");
                flint_printf("aa = ");
                fq_nmod_print_pretty(aa, ctx->fq_nmod_ctx);
                flint_printf("\n");
                flint_printf("a = ");
                fq_zech_print_pretty(a, ctx);
                flint_printf("\n");
                flint_printf("b = ");
                fq_zech_print_pretty(b, ctx);
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_clear(x);
            fq_nmod_clear(bb, ctx->fq_nmod_ctx);
            fq_nmod_clear(aa, ctx->fq_nmod_ctx);
        }

        fq_zech_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

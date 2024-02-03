/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"

TEST_FUNCTION_START(fq_zech_mul_fmpz, state)
{
    slong ix, jx;
    int result;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fq_zech_ctx_t ctx;

        fq_zech_ctx_init_randtest(ctx, state, 1);

        for (jx = 0; jx < 10; jx++)
        {
            fmpz_t x, p;
            fq_nmod_t aa, bb;
            fq_zech_t a, b, c;

            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_nmod_init(bb, ctx->fq_nmod_ctx);

            fmpz_init(x);
            fmpz_init_set_ui(p, fq_zech_ctx_prime(ctx));
            fmpz_randtest_mod_signed(x, state, p);
            fmpz_clear(p);

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

        for (jx = 0; jx < 10; jx++)
        {
            fmpz_t x, p;
            fq_nmod_t aa, bb;
            fq_zech_t a, b;

            fq_nmod_init(aa, ctx->fq_nmod_ctx);
            fq_nmod_init(bb, ctx->fq_nmod_ctx);

            fmpz_init(x);
            fmpz_init_set_ui(p, fq_zech_ctx_prime(ctx));
            fmpz_randtest_mod_signed(x, state, p);
            fmpz_clear(p);

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

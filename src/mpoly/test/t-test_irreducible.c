/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(mpoly_test_irreducible, state)
{
    slong i;

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;

        fmpz_mpoly_ctx_init(ctx, 8, ORD_LEX);
        fmpz_mpoly_init(f, ctx);

        fmpz_mpoly_set_str_pretty(f,
                "x1^809*x2^75*x3^384*x4^324*x5^74*x6^788*x7^83*x8^414+"
                "x1^805*x2^343*x3^595*x4^246*x5^32*x6^90*x7^473*x8^591+"
                "x1^718*x2^108*x3^680*x4^368*x5^358*x8^276+"
                "x1^683*x2^533*x4^649*x5^619*x6^136*x7^223*x8^610+"
                "x2^617*x3^777*x4^799*x5^443*x6^545*x7^166*x8^216+"
                "x1^485*x2^646*x3^424*x4^265*x5^416*x6^400*x7^278+"
                "x1^336*x2^149*x3^361*x4^691*x5^629*x6^282*x7^530*x8^259+"
                "x1^266*x3^258*x5^422*x6^637*x7^244*x8^236+"
                "x1^74*x2^812*x3^162*x4^417*x5^71*x6^188*x7^258*x8^637+"
                "x1^37*x2^604*x3^94*x4^474*x6^853*x7^521*x8^250", NULL, ctx);

        if (!mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo))
        {
            flint_printf("FAIL: check 8 variable example\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        const char * vars[] = {"x", "y"};

        fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpz_mpoly_init(f, ctx);

        fmpz_mpoly_set_str_pretty(f, "y^2639+x^4432*y^2436+x^400*y^1827+x^1300", vars, ctx);
        if (!mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo))
        {
            flint_printf("FAIL: check 2 variable example 1\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_set_str_pretty(f, "y^5481+x^2*y^5477+x^4167*y^4366+x^2700", vars, ctx);
        if (!mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo))
        {
            flint_printf("FAIL: check 2 variable example 1\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < 500*flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len1, len2;
        flint_bitcnt_t bits1, bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 10) + 2;
        len2 = n_randint(state, 10) + 2;
        bits1 = n_randint(state, 9) + 1;
        bits2 = n_randint(state, 9) + 1;

        fmpz_mpoly_randtest_bits(h, state, len1, 10, bits1, ctx);
        fmpz_mpoly_randtest_bits(g, state, len2, 10, bits2, ctx);
        fmpz_mpoly_mul(f, h, g, ctx);

        if (mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo) &&
            h->length > 1 && g->length > 1)
        {
            flint_printf("FAIL: check reducible input\n");
            flint_printf("f: "); fmpz_mpoly_print_pretty(f, NULL, ctx); flint_printf("\n");
            flint_printf("g: "); fmpz_mpoly_print_pretty(g, NULL, ctx); flint_printf("\n");
            flint_printf("h: "); fmpz_mpoly_print_pretty(h, NULL, ctx); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

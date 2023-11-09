/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly_factor.h"

/* Defined in t-gcd_brown.c, t-gcd_brown_threaded.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#define compute_gcd compute_gcd_zippel
int compute_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t state;
    int success = 0;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    fmpz_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    slong * perm;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }

        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 2)
    {
        return fmpz_mpoly_gcd(G, A, B, ctx);
    }

    perm = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);
    shift = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);
    stride = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    flint_randinit(state);

    wbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(lctx, ctx->minfo->nvars, ORD_LEX);

    fmpz_mpoly_init3(Al, 0, wbits, lctx);
    fmpz_mpoly_init3(Bl, 0, wbits, lctx);
    fmpz_mpoly_init3(Gl, 0, wbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, wbits, lctx);
    fmpz_mpoly_init3(Ac, 0, wbits, lctx);
    fmpz_mpoly_init3(Bc, 0, wbits, lctx);
    fmpz_mpoly_init3(Gc, 0, wbits, lctx);

    fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
    fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);

    success = fmpz_mpolyl_content(Ac, Al, 1, lctx) &&
              fmpz_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    if (!fmpz_mpoly_divides(Al, Al, Ac, lctx) ||
        !fmpz_mpoly_divides(Bl, Bl, Bc, lctx))
    {
        flint_printf("FAIL: Check content divides\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mpoly_repack_bits_inplace(Al, wbits, lctx);
    fmpz_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    success = fmpz_mpolyl_gcd_zippel(Gl, Abarl, Bbarl, Al, Bl, lctx, state);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    fmpz_mpoly_mul(Gl, Gl, Gc, lctx);
    fmpz_mpoly_repack_bits_inplace(Gl, wbits, lctx);

    fmpz_mpoly_from_mpolyl_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                                Gl, lctx, perm, shift, stride);
    success = 1;

    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

cleanup:

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_clear(Ac, lctx);
    fmpz_mpoly_clear(Bc, lctx);
    fmpz_mpoly_clear(Gc, lctx);

    fmpz_mpoly_ctx_clear(lctx);

    flint_randclear(state);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}

/* Defined in t-gcd_brown.c, t-gcd_brown_threaded.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#ifndef gcd_check
#define gcd_check gcd_check
void gcd_check(
    fmpz_mpoly_t g,
    fmpz_mpoly_t a,
    fmpz_mpoly_t b,
    const fmpz_mpoly_t gdiv,
    fmpz_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name,
    int compute_gcd_fun(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t))
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = compute_gcd_fun(g, a, b, ctx);

    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpz_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        flint_printf("FAIL: Check gcd has positive lc\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = compute_gcd_fun(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mpoly_clear(ca, ctx);
    fmpz_mpoly_clear(cb, ctx);
    fmpz_mpoly_clear(cg, ctx);
}
#endif

TEST_FUNCTION_START(fmpz_mpoly_factor_gcd_zippel, state)
{
    slong i, j, tmul = 25;

    /* examples from Zippel's 1979 paper */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t r, d, f, g;
        const char* vars[] =
            {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};

        const char * example[][3] =
        {{
            "x1^2 + x1 + 3",
            "2*x1^2 + 2*x1 + 1",
            "x1^2 + 2*x1 + 2"
        }, {
            "2*x1^2*x2^2 + x1*x2 + 2*x1",
            "x2^2 + 2*x1^2*x2 + x1^2 + 1",
            "x1^2*x2^2 + x1^2*x2 + x1*x2 + x1^2 + x1"
        }, {
            "x2^2*x3^2 + x2^2*x3 + 2*x1^2*x2*x3 + x1*x3",
            "x3^2 + x2^2*x3 + x1^2*x2*x3 + x1*x3 + x1^2*x2^2",
            "x2*x3 + 2*x1*x3 + x3 + x1"
        }, {
            "x1^2*x4^2 + x2^2*x3*x4 + x1^2*x2*x4 + x2*x4 + x1^2*x2*x3",
            "x1*x2*x3^2*x4^2 + x1*x3^2*x4^2 + x1*x4^2 + x4^2 + x1*x3*x4",
            "x1*x3^2*x4^2 + x3^2*x4^2 + x4^2 + x1*x2^2*x3*x4 + x1*x2^2"
        }, {
            "x1^3*x2^2*x3^2*x4*x5^2 + x1*x2^2*x5^2 + x1^3*x3*x4^2*x5"
                                    " + x1^3*x2*x3^2*x4*x5 + x1^2*x2*x3^2*x4^2"
            ,
            "x1*x2^2*x5^2 + x1*x2*x3^2*x4*x5 + x1*x2*x3^2*x4^2"
                                                          " + x1*x2^2*x4^2 + 1"
            ,
            "x1*x3^2*x4*x5^2 + x2*x5^2 + x1*x2*x4*x5 + x2*x5 + x1*x2*x3*x4^2"
        }, {
            "x1*x2*x4^2*x5^2*x6^2 + x1*x2^2*x3^2*x4*x5^2*x6^2  + x1^2*x3*x6^2"
                                   " + x1^2*x2*x3^2*x4*x5^2*x6 + x1^2*x3*x5*x6"
            ,
            "x1^2*x2*x4*x5^2*x6^2 + x1*x3*x5^2*x6^2 + x1*x2^2*x6^2"
                                      " + x1^2*x2^2*x3^2*x5*x6 + x1*x3^2*x4*x5"
            ,
            "x2^2*x3^2*x4*x5^2*x6 + x1*x4^2*x5*x6 + x2^2*x3^2*x4*x5*x6"
                                         " + x1*x2^2*x3*x4^2*x6 + x1^2*x3*x5^2"
        }, {
            "x1*x2^2*x4^2*x6^2*x7^2 + x1^2*x3*x4*x6^2*x7^2 + x3^2*x4^2*x7^2"
                                              " + x1^2*x2*x4^2*x6 + x3*x4*x5^2"
            ,
            "x1^2*x2*x4^2*x5*x6^2*x7^2 + x1*x2*x3*x6*x7 + x3*x4^2*x5^2*x7"
                                   " + x1*x4^2*x5^2*x7 + x1^2*x2*x3*x4^2+x5*x6"
            ,
            "x1*x3*x5*x6^2*x7^2 + x2^2*x3^2*x4^2*x5*x6*x7^2 + x4*x6*x7^2"
                                   " + x1^2*x2*x3*x5*x6*x7 + x1^2*x3^2*x4*x5^2"
        }, {
            "x2^2*x4*x5*x6*x7*x8^2 + x1^2*x2*x3^2*x4^2*x6^2*x7^2*x8"
                   " + x1^2*x3*x4^2*x6^2*x7^2 + x1^2*x2^2*x3^2*x4*x5^2*x6*x7^2"
                                                                " + x2^2*x4*x6"
            ,
            "x1^2*x2^2*x3*x4^2*x5*x6^2*x8^2 + x2*x5*x6^2*x8^2"
              " + x1^2*x2^2*x3^2*x4^2*x6^2*x7^2*x8 + x1^2*x3^2*x4*x5^2*x7^2*x8"
                                                      " + x1*x2^2*x3^2*x5^2*x7"
            ,
            "x1*x4^2*x5*x6*x7*x8^2 + x1*x2^2*x4^2*x5^2*x6^2*x8"
                       " + x1^2*x2*x3*x4^2*x6^2*x8 + x1^2*x2^2*x3^2*x4*x5^2*x8"
                                                           " + x1*x2*x4^2*x5^2"
        }, {
            "x1^2*x3^3*x4*x6*x8*x9^2 + x1*x2*x3*x4^2*x5^2*x8*x9"
                      " + x2*x3*x4*x5^2*x8*x9 + x1*x3^3*x4^2*x5^2*x6^2*x7*x8^2"
                                                  " + x2*x3*x4*x5^2*x6*x7*x8^2"
            ,
            "x1^2*x2^2*x3*x7^2*x8*x9 + x2^2*x9 + x1^2*x3*x4^2*x5^2*x6*x7^2"
                                            " + x4^2*x5^2*x7^2 + x3*x4^2*x6*x7"
            ,
            "x1^2*x2*x4*x5*x6*x7^2*x8^2*x9^2 + x1^2*x2*x3*x5*x6^2*x7^2*x8*x9^2"
                              " + x1^2*x3*x4*x6*x7^2*x8*x9 + x1^2*x2^2*x6*x8^2"
                                                        " + x2^2*x4*x5*x6^2*x7"
        }, {
            "x1*x2^2*x4^2*x8*x9^2*x10^2 + x2^2*x4*x5^2*x6*x7*x9*x10^2"
                        " + x1^2*x2*x3*x5^2*x7^2*x9^2 + x1*x3^2*x4^2*x7^2*x9^2"
                                                      " + x1^2*x3*x4*x7^2*x8^2"
            ,
            "x1*x2*x3^2*x4*x6*x7*x8*x9^2*x10^2 + x2^2*x3^2*x4^2*x6^2*x9*x10^2"
                                    " + x1*x2^2*x3^2*x4*x5*x6*x7*x8^2*x9^2*x10"
               " + x1^2*x2*x4^2*x5^2*x8^2*x9^2*x10 + x3*x4^2*x5*x6*x7^2*x9*x10"
            ,
            "x1*x2^2*x3^2*x5^2*x6^2*x7*x8*x9^2*x10^2 + x3*x8*x9^2*x10^2"
                  " + x1*x2^2*x3*x4*x5^2*x6^2*x8^2*x9*x10 + x1*x3*x6*x7*x8*x10"
                                                    " + x4^2*x5^2*x6^2*x7*x9^2"
        }};

        for (i = 1; i <= 10; i++)
        {
            fmpz_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX);
            fmpz_mpoly_init(r, ctx);
            fmpz_mpoly_init(d, ctx);
            fmpz_mpoly_init(f, ctx);
            fmpz_mpoly_init(g, ctx);
            fmpz_mpoly_set_str_pretty(d, example[i - 1][0], vars, ctx);
            fmpz_mpoly_set_str_pretty(f, example[i - 1][1], vars, ctx);
            fmpz_mpoly_set_str_pretty(g, example[i - 1][2], vars, ctx);
            fmpz_mpoly_mul(f, f, d, ctx);
            fmpz_mpoly_mul(g, g, d, ctx);

            fmpz_mpoly_randtest_bits(r, state, 10, 100, FLINT_BITS, ctx);
            gcd_check(r, f, g, d, ctx, -2, i, "example", compute_gcd);

            fmpz_mpoly_clear(r, ctx);
            fmpz_mpoly_clear(d, ctx);
            fmpz_mpoly_clear(f, ctx);
            fmpz_mpoly_clear(g, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }
    }

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b, c;
        const char * vars[] = {"t" ,"z", "y", "x"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(c, ctx);

        fmpz_mpoly_set_str_pretty(c, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_mul(a, a, c, ctx);
        fmpz_mpoly_mul(b, b, c, ctx);

        gcd_check(g, a, b, c, ctx, -1, 1, "example", compute_gcd);

        /*
            Main variable is t and c has content wrt t mod first chosen prime.
            However, this content doesn't depend on z and must be caught before
            it gets down to bivariate. Buggy code hangs on this example.
        */
        fmpz_mpoly_set_str_pretty(c, "t*(y+x)+9223372036854775837", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "1+x^2+y^2+z^2+t^2", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1+x^3+y^3+z^3+t^3", vars, ctx);
        fmpz_mpoly_mul(a, a, c, ctx);
        fmpz_mpoly_mul(b, b, c, ctx);

        gcd_check(g, a, b, c, ctx, -1, 1, "bug check", compute_gcd);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        if (ctx->minfo->nvars < 1)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 25) + 1;
        len1 = n_randint(state, 25) + 1;
        len2 = n_randint(state, 25) + 1;

        degbound = 2 + 100/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 300);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            if (fmpz_mpoly_is_zero(t, ctx))
                fmpz_mpoly_one(t, ctx);

            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse", compute_gcd);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef compute_gcd

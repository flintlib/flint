/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly_factor.h"

/* Defined in t-gcd_subresultant.c and t-gcd_zippel.c */
#define compute_gcd compute_gcd_zippel
int compute_gcd(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_bitcnt_t wbits;
    flint_rand_t state;
    int success = 0;
    nmod_mpoly_ctx_t lctx;
    nmod_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    nmod_mpoly_t Ac, Bc, Gc;
    ulong * shift, * stride;
    slong * perm;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
            nmod_mpoly_zero(G, ctx);
        else
            nmod_mpoly_make_monic(G, B, ctx);

        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
    {
        return 0;
    }

    if (ctx->minfo->nvars < 2)
    {
        return nmod_mpoly_gcd(G, A, B, ctx);
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

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(!nmod_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!nmod_mpoly_is_zero(B, ctx));
    FLINT_ASSERT(ctx->minfo->nvars >= 2);

    flint_randinit(state);

    wbits = FLINT_MAX(A->bits, B->bits);

    nmod_mpoly_ctx_init(lctx, ctx->minfo->nvars, ORD_LEX, nmod_mpoly_ctx_modulus(ctx));

    nmod_mpoly_init3(Al, 0, wbits, lctx);
    nmod_mpoly_init3(Bl, 0, wbits, lctx);
    nmod_mpoly_init3(Gl, 0, wbits, lctx);
    nmod_mpoly_init3(Abarl, 0, wbits, lctx);
    nmod_mpoly_init3(Bbarl, 0, wbits, lctx);
    nmod_mpoly_init3(Ac, 0, wbits, lctx);
    nmod_mpoly_init3(Bc, 0, wbits, lctx);
    nmod_mpoly_init3(Gc, 0, wbits, lctx);

    nmod_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
    nmod_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);

    success = nmod_mpolyl_content(Ac, Al, 1, lctx) &&
              nmod_mpolyl_content(Bc, Bl, 1, lctx);
    if (!success)
        goto cleanup;

    if (!nmod_mpoly_divides(Al, Al, Ac, lctx) ||
        !nmod_mpoly_divides(Bl, Bl, Bc, lctx))
    {
        flint_printf("FAIL: Check content divides\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_repack_bits_inplace(Al, wbits, lctx);
    nmod_mpoly_repack_bits_inplace(Bl, wbits, lctx);

    success = nmod_mpolyl_gcdp_zippel_smprime(Gl, Abarl, Bbarl, Al, Bl,
                                           ctx->minfo->nvars - 1, lctx, state);
    if (!success)
        goto cleanup;

    success = nmod_mpoly_gcd(Gc, Ac, Bc, lctx);
    if (!success)
        goto cleanup;

    nmod_mpoly_mul(Gl, Gl, Gc, lctx);

    nmod_mpoly_from_mpolyl_perm_inflate(G, FLINT_MIN(A->bits, B->bits), ctx,
                                                Gl, lctx, perm, shift, stride);
    success = 1;

    nmod_mpoly_make_monic(G, G, ctx);

cleanup:

    nmod_mpoly_clear(Al, lctx);
    nmod_mpoly_clear(Bl, lctx);
    nmod_mpoly_clear(Gl, lctx);
    nmod_mpoly_clear(Abarl, lctx);
    nmod_mpoly_clear(Bbarl, lctx);
    nmod_mpoly_clear(Ac, lctx);
    nmod_mpoly_clear(Bc, lctx);
    nmod_mpoly_clear(Gc, lctx);

    nmod_mpoly_ctx_clear(lctx);

    flint_randclear(state);

    flint_free(perm);
    flint_free(shift);
    flint_free(stride);

    return success;
}

/* Defined in t-gcd_subresultant.c and t-gcd_zippel.c */
#define gcd_check gcd_check_gcd_zippel
void gcd_check(
    nmod_mpoly_t g,
    nmod_mpoly_t a,
    nmod_mpoly_t b,
    const nmod_mpoly_t gdiv,
    nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    nmod_mpoly_t ca, cb, cg;

    nmod_mpoly_init(ca, ctx);
    nmod_mpoly_init(cb, ctx);
    nmod_mpoly_init(cg, ctx);

    res = compute_gcd(g, a, b, ctx);

    nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        if (FLINT_BIT_COUNT(nmod_mpoly_ctx_modulus(ctx)) < 10)
            goto cleanup;

        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_is_zero(gdiv, ctx))
    {
        if (!nmod_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (nmod_mpoly_is_zero(g, ctx))
    {
        if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (g->coeffs[0] != 1)
    {
        flint_printf("FAIL: Check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = nmod_mpoly_divides(ca, a, g, ctx) &&
          nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = compute_gcd(cg, ca, cb, ctx);
    nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        if (FLINT_BIT_COUNT(nmod_mpoly_ctx_modulus(ctx)) < 10)
            goto cleanup;

        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_factor_gcd_zippel, state)
{
    slong i, j, tmul = 20;

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t r, a, b, g;

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 1009);

        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(a, "x1+x2+x3+1", NULL, ctx);
        nmod_mpoly_set_str_pretty(b, "x1+x2+x3+2", NULL, ctx);
        nmod_mpoly_set_str_pretty(g, "x1^3*(x3+1)*x4 + x1^2*(x3+1)*(x4+1) + x1*(x3+1)*(x3+2)*x2 + (x3^2+x3+2)*x4 + 1", NULL, ctx);

        nmod_mpoly_mul(a, a, g, ctx);
        nmod_mpoly_mul(b, b, g, ctx);

        gcd_check(r, a, b, g, ctx, -2, 0, "example");

        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b, t;
        const char * vars[] = {"t" ,"z", "y", "x"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 1000003);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(t, ctx);

        nmod_mpoly_set_str_pretty(t, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        nmod_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        nmod_mpoly_mul(a, a, t, ctx);
        nmod_mpoly_mul(b, b, t, ctx);

        gcd_check(g, a, b, t, ctx, -2, 1, "example");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t r, d, f, g;
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
            nmod_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX, 1009);
            nmod_mpoly_init(r, ctx);
            nmod_mpoly_init(d, ctx);
            nmod_mpoly_init(f, ctx);
            nmod_mpoly_init(g, ctx);
            nmod_mpoly_set_str_pretty(d, example[i - 1][0], vars, ctx);
            nmod_mpoly_set_str_pretty(f, example[i - 1][1], vars, ctx);
            nmod_mpoly_set_str_pretty(g, example[i - 1][2], vars, ctx);
            nmod_mpoly_mul(f, f, d, ctx);
            nmod_mpoly_mul(g, g, d, ctx);

            nmod_mpoly_randtest_bits(r, state, 10, FLINT_BITS, ctx);
            gcd_check(r, f, g, d, ctx, -1, i, "example");

            nmod_mpoly_clear(r, ctx);
            nmod_mpoly_clear(d, ctx);
            nmod_mpoly_clear(f, ctx);
            nmod_mpoly_clear(g, ctx);
            nmod_mpoly_ctx_clear(ctx);
        }
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t p;

        p = n_randint(state, (i % 2 == 0) ? 10 : FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, p);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 30) + 1;
        len1 = n_randint(state, 30) + 1;
        len2 = n_randint(state, 30) + 1;

        degbound = 2 + 150/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);

            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);
            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef compute_gcd
#undef gcd_check

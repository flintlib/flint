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
#include "thread_support.h"
#include "fmpz_mpoly_factor.h"

typedef struct
{
    fmpz_mpoly_struct * Pl;
    const fmpz_mpoly_ctx_struct * lctx;
    const fmpz_mpoly_struct * P;
    const fmpz_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift;
    const ulong * stride;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertl_arg_struct;

typedef _convertl_arg_struct _convertl_arg_t[1];

static void _worker_convertu(void * varg)
{
    _convertl_arg_struct * arg = (_convertl_arg_struct *) varg;

    fmpz_mpoly_to_mpolyl_perm_deflate(arg->Pl, arg->lctx, arg->P, arg->ctx,
                                         arg->perm, arg->shift, arg->stride);
}

/* Defined in t-gcd_brown.c, t-gcd_brown_threaded.c, t-gcd_subresultant.c,
 * t-gcd_zippel.c, t-gcd_zippel2.c */
#define compute_gcd compute_gcd_brown_threaded
int compute_gcd(
    fmpz_mpoly_t G,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong * perm;
    ulong * shift, * stride;
    slong i;
    flint_bitcnt_t ABbits;
    fmpz_mpoly_ctx_t lctx;
    fmpz_mpoly_t Al, Bl, Gl, Abarl, Bbarl;
    thread_pool_handle * handles;
    slong num_handles;
    slong thread_limit = FLINT_MIN(A->length, B->length)/16;

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

    perm = (slong *) flint_malloc(ctx->minfo->nvars*sizeof(slong));
    shift = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    stride = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        perm[i] = i;
        shift[i] = 0;
        stride[i] = 1;
    }

    ABbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(lctx, ctx->minfo->nvars, ORD_LEX);
    fmpz_mpoly_init3(Al, 0, ABbits, lctx);
    fmpz_mpoly_init3(Bl, 0, ABbits, lctx);
    fmpz_mpoly_init3(Gl, 0, ABbits, lctx);
    fmpz_mpoly_init3(Abarl, 0, ABbits, lctx);
    fmpz_mpoly_init3(Bbarl, 0, ABbits, lctx);

    num_handles = flint_request_threads(&handles, thread_limit);

    /* convert inputs */
    if (num_handles > 0)
    {
        slong m = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertl_arg_t arg;

        FLINT_ASSERT(m >= 0);
        FLINT_ASSERT(m < num_handles);

        arg->Pl = Bl;
        arg->lctx = lctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = perm;
        arg->shift = shift;
        arg->stride = stride;
        arg->handles = handles + (m + 1);
        arg->num_handles = num_handles - (m + 1);

        thread_pool_wake(global_thread_pool, handles[m], 0, _worker_convertu, arg);

        fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);

        thread_pool_wait(global_thread_pool, handles[m]);
    }
    else
    {
        fmpz_mpoly_to_mpolyl_perm_deflate(Al, lctx, A, ctx, perm, shift, stride);
        fmpz_mpoly_to_mpolyl_perm_deflate(Bl, lctx, B, ctx, perm, shift, stride);
    }

    /* calculate gcd */
    success = fmpz_mpolyl_gcd_brown_threaded_pool(Gl, Abarl, Bbarl, Al, Bl,
                                             lctx, NULL, handles, num_handles);

    flint_give_back_threads(handles, num_handles);

    if (success)
    {
        fmpz_mpoly_from_mpolyl_perm_inflate(G, ABbits, ctx, Gl, lctx,
                                                          perm, shift, stride);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
    }

    fmpz_mpoly_clear(Al, lctx);
    fmpz_mpoly_clear(Bl, lctx);
    fmpz_mpoly_clear(Gl, lctx);
    fmpz_mpoly_clear(Abarl, lctx);
    fmpz_mpoly_clear(Bbarl, lctx);
    fmpz_mpoly_ctx_clear(lctx);

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

TEST_FUNCTION_START(fmpz_mpoly_factor_gcd_brown_threaded, state)
{
    slong i, j, tmul = 12, max_threads = 5;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);
        if (ctx->minfo->nvars < 1)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 50) + 1;
        len1 = n_randint(state, 80);
        len2 = n_randint(state, 80);

        degbound = 1 + 50/ctx->minfo->nvars/ctx->minfo->nvars;

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
            fmpz_mpoly_scalar_mul_ui(a, a, n_randint(state, 10) + 1, ctx);
            fmpz_mpoly_scalar_mul_ui(b, b, n_randint(state, 10) + 1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "random dense", compute_gcd);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
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

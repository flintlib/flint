/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arf.h"
#include "gr.h"
#include "gr_vec.h"
#include "nfloat.h"

/* Checks nfloat_div, nfloat_inv, nfloat_div_ui, nfloat_div_si, nfloat_sqrt,
   nfloat_rsqrt and the vector division functions against the correctly
   rounded arf functions: with directed rounding the result must be the
   correctly rounded floor or ceiling, otherwise within 1 ulp of the exact
   value. Inputs include exact and nearly exact quotients and roots, which
   exercise the certification fallbacks. */

#define OP_DIV 0
#define OP_INV 1
#define OP_SQRT 2
#define OP_RSQRT 3
#define OP_DIV_UI 4
#define OP_DIV_SI 5

static const char * op_names[] = { "div", "inv", "sqrt", "rsqrt", "div_ui", "div_si" };

static int
nfloat_identical(nfloat_srcptr a, nfloat_srcptr b, gr_ctx_t ctx)
{
    if (NFLOAT_EXP(a) != NFLOAT_EXP(b))
        return 0;
    if (NFLOAT_IS_SPECIAL(a))
        return 1;
    return NFLOAT_SGNBIT(a) == NFLOAT_SGNBIT(b) &&
        flint_mpn_equal_p(NFLOAT_D(a), NFLOAT_D(b), NFLOAT_CTX_NLIMBS(ctx));
}

/* Random normal nfloat; mode selects the distribution. */
static void
nfloat_randtest_normal(nfloat_ptr x, flint_rand_t state, slong emax, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);

    do {
        if (n_randint(state, 4) == 0)
        {
            GR_IGNORE(nfloat_randtest(x, state, ctx));
        }
        else
        {
            for (i = 0; i < n; i++)
                NFLOAT_D(x)[i] = n_randtest(state);
            /* occasionally leave trailing zero limbs, or a power of two */
            if (n_randint(state, 8) == 0)
                flint_mpn_zero(NFLOAT_D(x), n_randint(state, n));
            if (n_randint(state, 32) == 0)
                flint_mpn_zero(NFLOAT_D(x), n);
            NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
            NFLOAT_EXP(x) = (slong) n_randint(state, 2 * emax + 1) - emax;
            NFLOAT_SGNBIT(x) = n_randint(state, 2);
        }
    }
    while (NFLOAT_IS_SPECIAL(x));
}

static void
print_limbs(const char * name, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    flint_printf("%s: exp = %wd, sgn = %wu, limbs (high to low) =", name, NFLOAT_EXP(x), NFLOAT_SGNBIT(x));
    for (i = n - 1; i >= 0; i--)
        flint_printf(" %016wx", NFLOAT_D(x)[i]);
    flint_printf("\n");
}

/* res = op(x, y, c) rounded to prec bits in the direction rnd (correctly
   rounded); y is used only by OP_DIV */
static void
ref_op(arf_t res, int op, nfloat_srcptr x, nfloat_srcptr y, ulong c, slong prec, arf_rnd_t rnd, gr_ctx_t ctx)
{
    arf_t xf, yf;

    arf_init(xf);
    arf_init(yf);
    GR_MUST_SUCCEED(nfloat_get_arf(xf, x, ctx));
    if (op == OP_DIV)
        GR_MUST_SUCCEED(nfloat_get_arf(yf, y, ctx));

    switch (op)
    {
        case OP_DIV: arf_div(res, xf, yf, prec, rnd); break;
        case OP_INV: arf_ui_div(res, 1, xf, prec, rnd); break;
        case OP_SQRT: arf_sqrt(res, xf, prec, rnd); break;
        case OP_RSQRT: arf_rsqrt(res, xf, prec, rnd); break;
        case OP_DIV_UI: arf_div_ui(res, xf, c, prec, rnd); break;
        default: arf_div_si(res, xf, (slong) c, prec, rnd); break;
    }

    arf_clear(xf);
    arf_clear(yf);
}

/* Compares z (computed with status) against op(x, y, c). */
static void
check_result(int op, nfloat_srcptr z, int status, gr_ctx_t ctx,
    nfloat_srcptr x, nfloat_srcptr y, ulong c)
{
    slong prec = NFLOAT_CTX_PREC(ctx);
    int flags = NFLOAT_CTX_FLAGS(ctx);
    arf_t exact, zf, ref, err, ulp;
    int fail = 0;

    arf_init(exact);
    arf_init(zf);
    arf_init(ref);
    arf_init(err);
    arf_init(ulp);

    /* the exact value truncated to 4 prec + 128 bits */
    ref_op(exact, op, x, y, c, 4 * prec + 128, ARF_RND_DOWN, ctx);

    if (status != GR_SUCCESS || NFLOAT_IS_SPECIAL(z))
    {
        /* only acceptable for results outside the exponent range,
           2^(e-1) <= |exact| < 2^e with e >= NFLOAT_MAX_EXP - 1 or
           e <= NFLOAT_MIN_EXP + 1 */
        int out_of_range = arf_is_special(exact) ||
            arf_cmpabs_2exp_si(exact, NFLOAT_MAX_EXP - 2) >= 0 ||
            arf_cmpabs_2exp_si(exact, NFLOAT_MIN_EXP + 1) < 0;

        if (!out_of_range)
        {
            flint_printf("FAIL (%s): unexpected status %d or special value\n", op_names[op], status);
            goto failure;
        }
        goto cleanup;
    }

    if (!LIMB_MSB_IS_SET(NFLOAT_D(z)[NFLOAT_CTX_NLIMBS(ctx) - 1]))
    {
        flint_printf("FAIL (%s): not normalized\n", op_names[op]);
        fail = 1;
    }

    GR_MUST_SUCCEED(nfloat_get_arf(zf, z, ctx));

    if (flags & NFLOAT_RND_FLOOR)
    {
        ref_op(ref, op, x, y, c, prec, ARF_RND_FLOOR, ctx);
        if (!arf_equal(zf, ref))
        {
            flint_printf("FAIL (%s): floor\n", op_names[op]);
            fail = 1;
        }
    }
    else if (flags & NFLOAT_RND_CEIL)
    {
        ref_op(ref, op, x, y, c, prec, ARF_RND_CEIL, ctx);
        if (!arf_equal(zf, ref))
        {
            flint_printf("FAIL (%s): ceil\n", op_names[op]);
            fail = 1;
        }
    }
    else
    {
        /* |z - exact| <= ulp(z) (1 + 2^-30), with exact known to
           4 prec + 128 bits */
        arf_sub(err, zf, exact, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(err, err);
        arf_one(ulp);
        arf_mul_2exp_si(ulp, ulp, NFLOAT_EXP(z) - prec);
        if (arf_cmp(err, ulp) > 0)
        {
            arf_div(err, err, ulp, 53, ARF_RND_UP);
            if (arf_cmp_d(err, 1.0 + ldexp(1.0, -30)) > 0)
            {
                flint_printf("FAIL (%s): error = %g ulp\n", op_names[op], arf_get_d(err, ARF_RND_UP));
                fail = 1;
            }
        }
    }

    if (!fail)
        goto cleanup;

failure:
    flint_printf("prec = %wd, flags = %d\n", prec, flags);
    flint_printf("x = %{gr}\n", x, ctx);
    if (op == OP_DIV)
        flint_printf("y = %{gr}\n", y, ctx);
    if (op == OP_DIV_UI)
        flint_printf("c = %wu\n", c);
    if (op == OP_DIV_SI)
        flint_printf("c = %wd\n", (slong) c);
    flint_printf("z = %{gr}\n", z, ctx);
    flint_printf("exact = "); arf_printd(exact, 40); flint_printf("\n");
    print_limbs("x", x, ctx);
    if (op == OP_DIV)
        print_limbs("y", y, ctx);
    print_limbs("z", z, ctx);
    flint_abort();

cleanup:
    arf_clear(exact);
    arf_clear(zf);
    arf_clear(ref);
    arf_clear(err);
    arf_clear(ulp);
}

static int
apply_op(int op, nfloat_ptr z, nfloat_srcptr x, nfloat_srcptr y, ulong c, gr_ctx_t ctx)
{
    switch (op)
    {
        case OP_DIV: return nfloat_div(z, x, y, ctx);
        case OP_INV: return nfloat_inv(z, x, ctx);
        case OP_SQRT: return nfloat_sqrt(z, x, ctx);
        case OP_RSQRT: return nfloat_rsqrt(z, x, ctx);
        case OP_DIV_UI: return nfloat_div_ui(z, x, c, ctx);
        default: return nfloat_div_si(z, x, (slong) c, ctx);
    }
}

TEST_FUNCTION_START(nfloat_div_sqrt, state)
{
    slong iter;

    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n, prec, emax;
        int flags, op, rnd, alias;
        ulong c = 0;
        nfloat_ptr x, y, z, w;
        int status;

        if (n_randint(state, 4) == 0)
            n = 1 + n_randint(state, NFLOAT_MAX_LIMBS);
        else
            n = 1 + n_randint(state, 10);

        prec = n * FLINT_BITS;
        rnd = n_randint(state, 3);
        flags = (rnd == 0) ? 0 : (rnd == 1) ? NFLOAT_RND_FLOOR : NFLOAT_RND_CEIL;
        if (n_randint(state, 4) == 0)
            flags |= NFLOAT_ALLOW_UNDERFLOW | NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN;

        nfloat_ctx_init(ctx, prec, flags);

        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        w = gr_heap_init(ctx);

        emax = (n_randint(state, 50) == 0) ? NFLOAT_MAX_EXP : 20;

        op = n_randint(state, 6);
        nfloat_randtest_normal(x, state, emax, ctx);
        nfloat_randtest_normal(y, state, emax, ctx);

        if (op == OP_SQRT || op == OP_RSQRT)
            NFLOAT_SGNBIT(x) = 0;

        if (op == OP_DIV_UI || op == OP_DIV_SI)
        {
            c = n_randtest_not_zero(state);
            if (op == OP_DIV_SI && c == (UWORD(1) << (FLINT_BITS - 1)))
                c = 3;
        }

        /* hard cases: nearly exact or exact results */
        switch (n_randint(state, 4))
        {
            case 0:
                if (op == OP_DIV)
                {
                    /* x = y q (rounded), q random */
                    nfloat_randtest_normal(w, state, 10, ctx);
                    if (n_randint(state, 2))
                    {
                        /* exact: q small integer, y with low bits zero */
                        GR_MUST_SUCCEED(nfloat_set_si(w, (slong) n_randint(state, 1000) - 500 + 1001 * (n_randint(state, 2) ? 1 : -1), ctx));
                        NFLOAT_D(y)[0] &= ~UWORD(0xfff);
                    }
                    status = nfloat_mul(x, y, w, ctx);
                    if (status != GR_SUCCESS)
                        nfloat_randtest_normal(x, state, 20, ctx);
                }
                else if (op == OP_SQRT)
                {
                    /* x = s^2, exact if s has few bits */
                    nfloat_randtest_normal(w, state, 20, ctx);
                    if (n_randint(state, 2))
                    {
                        /* keep prec / 2 bits, so that the square is exact */
                        slong k, drop = prec - prec / 2;
                        for (k = 0; k < drop / FLINT_BITS; k++)
                            NFLOAT_D(w)[k] = 0;
                        if (drop % FLINT_BITS)
                            NFLOAT_D(w)[k] &= ~((UWORD(1) << (drop % FLINT_BITS)) - 1);
                    }
                    GR_MUST_SUCCEED(nfloat_sqr(x, w, ctx));
                }
                else if (op == OP_RSQRT || op == OP_INV)
                {
                    /* x = 1 / s^2 resp. 1 / s, rounded */
                    nfloat_randtest_normal(w, state, 20, ctx);
                    NFLOAT_SGNBIT(w) = 0;
                    if (op == OP_RSQRT)
                        GR_MUST_SUCCEED(nfloat_sqr(w, w, ctx));
                    GR_IGNORE(nfloat_inv(x, w, ctx));
                    if (NFLOAT_IS_SPECIAL(x))
                        nfloat_randtest_normal(x, state, 20, ctx);
                    if (op == OP_RSQRT)
                        NFLOAT_SGNBIT(x) = 0;
                }
                else
                {
                    /* x = c q, q random */
                    nfloat_randtest_normal(w, state, 20, ctx);
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(nfloat_set_si(w, n_randint(state, 1000), ctx));
                    if (op == OP_DIV_UI)
                        status = gr_mul_ui(x, w, c, ctx);
                    else
                        status = gr_mul_si(x, w, (slong) c, ctx);
                    if (status != GR_SUCCESS || NFLOAT_IS_SPECIAL(x))
                        nfloat_randtest_normal(x, state, 20, ctx);
                }
                break;
            case 1:
                /* x == y */
                if (op == OP_DIV && n_randint(state, 4) == 0)
                    GR_MUST_SUCCEED(nfloat_set(y, x, ctx));
                break;
            default:
                break;
        }

        /* aliasing */
        alias = n_randint(state, 3);

        if (alias == 0 || op != OP_DIV)
        {
            if (alias == 0)
            {
                status = apply_op(op, z, x, y, c, ctx);
            }
            else
            {
                GR_MUST_SUCCEED(nfloat_set(z, x, ctx));
                status = apply_op(op, z, z, y, c, ctx);
            }
        }
        else if (alias == 1)
        {
            GR_MUST_SUCCEED(nfloat_set(z, x, ctx));
            status = apply_op(op, z, z, y, c, ctx);
        }
        else
        {
            GR_MUST_SUCCEED(nfloat_set(z, y, ctx));
            status = apply_op(op, z, x, z, c, ctx);
        }

        check_result(op, z, status, ctx, x, y, c);

        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_heap_clear(w, ctx);
        gr_ctx_clear(ctx);
    }

    /* vector functions */
    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n, prec, len, i, sz;
        int flags, rnd, which, inplace;
        gr_ptr X, Y, Z, W, c;
        ulong cu;
        int status;

        if (n_randint(state, 4) == 0)
            n = 1 + n_randint(state, NFLOAT_MAX_LIMBS);
        else
            n = 1 + n_randint(state, 10);

        prec = n * FLINT_BITS;
        rnd = n_randint(state, 3);
        flags = (rnd == 0) ? 0 : (rnd == 1) ? NFLOAT_RND_FLOOR : NFLOAT_RND_CEIL;
        nfloat_ctx_init(ctx, prec, flags);
        sz = ctx->sizeof_elem;

        len = n_randint(state, 10);
        X = gr_heap_init_vec(len, ctx);
        Y = gr_heap_init_vec(len, ctx);
        Z = gr_heap_init_vec(len, ctx);
        W = gr_heap_init_vec(len, ctx);
        c = gr_heap_init(ctx);

        for (i = 0; i < len; i++)
        {
            nfloat_randtest_normal(GR_ENTRY(X, i, sz), state, 20, ctx);
            nfloat_randtest_normal(GR_ENTRY(Y, i, sz), state, 20, ctx);
            if (n_randint(state, 10) == 0)
                GR_MUST_SUCCEED(nfloat_zero(GR_ENTRY(X, i, sz), ctx));
        }
        nfloat_randtest_normal(c, state, 20, ctx);
        if (n_randint(state, 8) == 0)
            flint_mpn_zero(NFLOAT_D(c), n - 1);

        /* nearly exact quotients for the scalar version */
        if (n_randint(state, 2))
        {
            for (i = 0; i < len; i++)
            {
                if (n_randint(state, 2))
                {
                    nn_ptr t = GR_ENTRY(W, i, sz);
                    nfloat_randtest_normal(t, state, 10, ctx);
                    if (n_randint(state, 2))
                        GR_MUST_SUCCEED(nfloat_set_si(t, (slong) n_randint(state, 100) - 50, ctx));
                    if (nfloat_mul(GR_ENTRY(X, i, sz), c, t, ctx) != GR_SUCCESS)
                        GR_MUST_SUCCEED(nfloat_zero(GR_ENTRY(X, i, sz), ctx));
                }
            }
        }

        cu = n_randtest_not_zero(state);
        which = n_randint(state, 4);
        inplace = n_randint(state, 2);

        if (inplace)
            GR_MUST_SUCCEED(_gr_vec_set(Z, X, len, ctx));

        if (which == 0)
            status = _gr_vec_div(Z, inplace ? Z : X, Y, len, ctx);
        else if (which == 1)
            status = _gr_vec_div_scalar(Z, inplace ? Z : X, len, c, ctx);
        else if (which == 2)
            status = _gr_vec_div_scalar_ui(Z, inplace ? Z : X, len, cu, ctx);
        else
            status = _gr_vec_div_scalar_si(Z, inplace ? Z : X, len, (slong) cu, ctx);

        for (i = 0; i < len; i++)
        {
            nn_ptr xi = GR_ENTRY(X, i, sz);
            nn_ptr yi = GR_ENTRY(Y, i, sz);
            nn_ptr zi = GR_ENTRY(Z, i, sz);
            int op;

            if (which == 0)
            {
                op = OP_DIV;
                /* identical to the scalar function */
                GR_MUST_SUCCEED(nfloat_div(GR_ENTRY(W, i, sz), xi, yi, ctx));
                if (!nfloat_identical(zi, GR_ENTRY(W, i, sz), ctx))
                {
                    flint_printf("FAIL: vec_div != div\n");
                    flint_abort();
                }
            }
            else if (which == 1)
            {
                op = OP_DIV;
                yi = c;
            }
            else
            {
                op = (which == 2) ? OP_DIV_UI : OP_DIV_SI;
            }

            if (NFLOAT_IS_ZERO(xi))
            {
                if (!NFLOAT_IS_ZERO(zi))
                {
                    flint_printf("FAIL: zero\n");
                    flint_abort();
                }
                continue;
            }

            check_result(op, zi, status, ctx, xi, yi, cu);
        }

        gr_heap_clear_vec(X, len, ctx);
        gr_heap_clear_vec(Y, len, ctx);
        gr_heap_clear_vec(Z, len, ctx);
        gr_heap_clear_vec(W, len, ctx);
        gr_heap_clear(c, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

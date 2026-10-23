/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "gr.h"
#include "gr_special.h"
#include "nfloat.h"

/* Checks the real elementary functions against ball evaluations in arb:
   in the default rounding mode the error must be at most about one ulp
   (relative to the function precision of the context), with directed
   rounding the result must be a valid bound within a few ulps. */

typedef struct
{
    const char * name;
    int method;
    int arity;
}
elem_func_t;

static const elem_func_t elem_funcs[] =
{
    { "exp", GR_METHOD_EXP, 1 },
    { "expm1", GR_METHOD_EXPM1, 1 },
    { "exp2", GR_METHOD_EXP2, 1 },
    { "log", GR_METHOD_LOG, 1 },
    { "log1p", GR_METHOD_LOG1P, 1 },
    { "log2", GR_METHOD_LOG2, 1 },
    { "sin", GR_METHOD_SIN, 1 },
    { "cos", GR_METHOD_COS, 1 },
    { "tan", GR_METHOD_TAN, 1 },
    { "sin_pi", GR_METHOD_SIN_PI, 1 },
    { "cos_pi", GR_METHOD_COS_PI, 1 },
    { "tan_pi", GR_METHOD_TAN_PI, 1 },
    { "atan", GR_METHOD_ATAN, 1 },
    { "atan2", GR_METHOD_ATAN2, 2 },
    { "sinh", GR_METHOD_SINH, 1 },
    { "cosh", GR_METHOD_COSH, 1 },
    { "tanh", GR_METHOD_TANH, 1 },
    { "pow", GR_METHOD_POW, 2 },
    { "exp10", GR_METHOD_EXP10, 1 },
    { "log10", GR_METHOD_LOG10, 1 },
    { "cot", GR_METHOD_COT, 1 },
    { "sec", GR_METHOD_SEC, 1 },
    { "csc", GR_METHOD_CSC, 1 },
    { "sinc", GR_METHOD_SINC, 1 },
    { "cot_pi", GR_METHOD_COT_PI, 1 },
    { "sec_pi", GR_METHOD_SEC_PI, 1 },
    { "csc_pi", GR_METHOD_CSC_PI, 1 },
    { "sinc_pi", GR_METHOD_SINC_PI, 1 },
    { "coth", GR_METHOD_COTH, 1 },
    { "sech", GR_METHOD_SECH, 1 },
    { "csch", GR_METHOD_CSCH, 1 },
    { "asin", GR_METHOD_ASIN, 1 },
    { "acos", GR_METHOD_ACOS, 1 },
    { "asin_pi", GR_METHOD_ASIN_PI, 1 },
    { "acos_pi", GR_METHOD_ACOS_PI, 1 },
    { "atan_pi", GR_METHOD_ATAN_PI, 1 },
    { "acot", GR_METHOD_ACOT, 1 },
    { "asec", GR_METHOD_ASEC, 1 },
    { "acsc", GR_METHOD_ACSC, 1 },
    { "acot_pi", GR_METHOD_ACOT_PI, 1 },
    { "asec_pi", GR_METHOD_ASEC_PI, 1 },
    { "acsc_pi", GR_METHOD_ACSC_PI, 1 },
    { "asinh", GR_METHOD_ASINH, 1 },
    { "acosh", GR_METHOD_ACOSH, 1 },
    { "atanh", GR_METHOD_ATANH, 1 },
    { "acoth", GR_METHOD_ACOTH, 1 },
    { "asech", GR_METHOD_ASECH, 1 },
    { "acsch", GR_METHOD_ACSCH, 1 },
    { "hypot", GR_METHOD_HYPOT, 2 },
};

#define NUM_FUNCS (sizeof(elem_funcs) / sizeof(elem_func_t))

static int
call_func(const elem_func_t * f, gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    if (f->arity == 1)
        return ((gr_method_unary_op) ctx->methods[f->method])(res, x, ctx);
    else
        return ((gr_method_binary_op) ctx->methods[f->method])(res, x, y, ctx);
}

/* random normal nfloat with interesting values */
static void
randtest_arg(nfloat_ptr x, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong r = n_randint(state, 16);

    for (i = 0; i < n; i++)
        NFLOAT_D(x)[i] = n_randtest(state);
    NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    NFLOAT_SGNBIT(x) = n_randint(state, 2);

    if (r == 0)
    {
        /* powers of two and small integers */
        GR_MUST_SUCCEED(nfloat_set_si(x, (slong) n_randint(state, 21) - 10, ctx));
        if (NFLOAT_IS_ZERO(x))
            GR_MUST_SUCCEED(nfloat_one(x, ctx));
        NFLOAT_EXP(x) += (slong) n_randint(state, 21) - 10;
    }
    else if (r == 1)
    {
        /* close to +-1 */
        slong k = n_randint(state, FLINT_BITS * n + 20);
        GR_MUST_SUCCEED(nfloat_one(x, ctx));
        if (n_randint(state, 2))
            GR_MUST_SUCCEED(nfloat_neg(x, x, ctx));
        {
            ulong t[NFLOAT_MAX_ALLOC];
            nfloat_set(t, x, ctx);
            NFLOAT_EXP(t) -= k;
            if (n_randint(state, 2))
                NFLOAT_SGNBIT(t) ^= 1;
            GR_IGNORE(nfloat_add(x, x, t, ctx));
        }
    }
    else if (r == 2)
    {
        /* huge or tiny */
        if (n_randint(state, 2))
            NFLOAT_EXP(x) = (slong) n_randint(state, 200) - 100;
        else if (n_randint(state, 4) == 0)
            /* the thresholds of exp, exp10, sech, csch, and huge trig
               arguments reduced by mp_real */
            NFLOAT_EXP(x) = n_randint(state, 2) ? FLINT_BITS - 10 + (slong) n_randint(state, 10)
                                                : 65536 + (slong) n_randint(state, 1000000);
        else if (n_randint(state, 2))
            NFLOAT_EXP(x) = (slong) n_randint(state, 2 * FLINT_BITS * n + 200) - FLINT_BITS * n - 100;
        else
            NFLOAT_EXP(x) = (n_randint(state, 2) ? 1 : -1) * (slong) n_randint(state, NFLOAT_MAX_EXP);
    }
    else if (r == 3)
    {
        /* few bits */
        for (i = 0; i < n - 1; i++)
            NFLOAT_D(x)[i] = 0;
        NFLOAT_D(x)[n - 1] &= ~((UWORD(1) << n_randint(state, FLINT_BITS)) - 1);
        NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
        NFLOAT_EXP(x) = (slong) n_randint(state, 40) - 20;
    }
    else if (r == 4)
    {
        /* close to a multiple of pi/2 (or pi/4) */
        ulong k = 1 + n_randtest(state) % (n_randint(state, 2) ? 10 : 1000000000);
        GR_MUST_SUCCEED(nfloat_pi(x, ctx));
        GR_IGNORE(gr_mul_ui(x, x, k, ctx));
        NFLOAT_EXP(x) -= 1 + n_randint(state, 2);
        if (n_randint(state, 2))
            NFLOAT_SGNBIT(x) = 1;
        if (n_randint(state, 2))
            NFLOAT_D(x)[0] ^= n_randint(state, 4);
    }
    else if (r == 5)
    {
        /* small, with the leading terms of series cancelling */
        NFLOAT_EXP(x) = -(slong) n_randint(state, FLINT_BITS / 2 * n + 20);
    }
    else
    {
        NFLOAT_EXP(x) = (slong) n_randint(state, 30) - 15;
    }
}

/* the reference value as an arb ball, at increasing precision until the
   radius is below 2^-30 ulp of the nfloat result (or the result is
   exact); returns 0 if arb fails */
static int
reference(arb_t ref, const elem_func_t * f, nfloat_srcptr x, nfloat_srcptr y, slong ulp_exp, gr_ctx_t ctx)
{
    gr_ctx_t actx;
    arb_t ax, ay;
    slong prec = NFLOAT_CTX_PREC(ctx), wp, maxwp;
    int status, ok = 0;

    arb_init(ax);
    arb_init(ay);
    GR_MUST_SUCCEED(nfloat_get_arf(arb_midref(ax), x, ctx));
    if (y != NULL)
        GR_MUST_SUCCEED(nfloat_get_arf(arb_midref(ay), y, ctx));

    /* periodic in x with period 2 (or 1): large x are even integers */
    if ((f->method == GR_METHOD_SIN_PI || f->method == GR_METHOD_COS_PI ||
        f->method == GR_METHOD_TAN_PI || f->method == GR_METHOD_COT_PI ||
        f->method == GR_METHOD_SEC_PI || f->method == GR_METHOD_CSC_PI) &&
        !NFLOAT_IS_SPECIAL(x) &&
        NFLOAT_EXP(x) > prec + 1)
        arb_zero(ax);

    /* (arb reduces trigonometric arguments up to 2^max(65536, 4 wp),
       as the slow path of nfloat does) */
    maxwp = 32 * prec + 5000;
    if (!NFLOAT_IS_SPECIAL(x) && NFLOAT_EXP(x) > 0 && NFLOAT_EXP(x) < (WORD(1) << 22))
        maxwp = FLINT_MAX(maxwp, NFLOAT_EXP(x) / 2);

    for (wp = 2 * prec + 64; wp < maxwp; wp *= 2)
    {
        gr_ctx_init_real_arb(actx, wp);
        if (f->method == GR_METHOD_HYPOT)
        {
            /* (no generic gr method) */
            arb_hypot(ref, ax, ay, wp);
            status = GR_SUCCESS;
        }
        else
            status = call_func(f, ref, ax, ay, actx);
        gr_ctx_clear(actx);

        if (status == GR_DOMAIN)
            break;

        /* (also a huge argument beyond arb's reduction at this precision) */
        if (status != GR_SUCCESS || !arb_is_finite(ref))
            continue;

        if (arb_is_exact(ref) || mag_cmp_2exp_si(arb_radref(ref), ulp_exp - 30) < 0)
        {
            ok = 1;
            break;
        }
    }

    arb_clear(ax);
    arb_clear(ay);
    return ok;
}

TEST_FUNCTION_START(nfloat_elem, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n, prec, fprec;
        int flags, rnd;
        const elem_func_t * f;
        nfloat_ptr x, y, z;
        arb_t ref, t;
        arf_t zf, ulp;
        int status, aliasing;
        slong ulp_exp;

        if (n_randint(state, 16) == 0)
            n = NFLOAT_MAX_LIMBS;   /* (the slow paths) */
        else if (n_randint(state, 4) == 0)
            n = 1 + n_randint(state, NFLOAT_MAX_LIMBS);
        else
            n = 1 + n_randint(state, 5);

        prec = n * FLINT_BITS;
        rnd = n_randint(state, 3);
        flags = (rnd == 0) ? 0 : (rnd == 1) ? NFLOAT_RND_FLOOR : NFLOAT_RND_CEIL;
        nfloat_ctx_init(ctx, prec, flags);

        fprec = prec;
        if (n_randint(state, 4) == 0)
        {
            fprec = 2 + n_randint(state, prec - 1);
            GR_MUST_SUCCEED(nfloat_ctx_set_func_prec(ctx, fprec));
        }

        f = elem_funcs + n_randint(state, NUM_FUNCS);

        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        arb_init(ref);
        arb_init(t);
        arf_init(zf);
        arf_init(ulp);

        randtest_arg(x, state, ctx);
        randtest_arg(y, state, ctx);

        aliasing = n_randint(state, 2);
        if (aliasing)
        {
            GR_MUST_SUCCEED(nfloat_set(z, x, ctx));
            status = call_func(f, z, z, y, ctx);
        }
        else
            status = call_func(f, z, x, y, ctx);

        if (status == GR_SUCCESS && !NFLOAT_IS_SPECIAL(z))
        {
            ulp_exp = NFLOAT_EXP(z) - prec;

            if (!reference(ref, f, x, (f->arity == 2) ? y : NULL, ulp_exp, ctx))
            {
                flint_printf("FAIL (%s): reference failed but nfloat succeeded\n", f->name);
                goto fail;
            }

            GR_MUST_SUCCEED(nfloat_get_arf(zf, z, ctx));

            /* t = z - ref */
            arb_sub_arf(t, ref, zf, ARF_PREC_EXACT);
            arb_neg(t, t);

            if (flags & NFLOAT_RND_FLOOR)
            {
                /* z <= ref, and within 3 ulps (of the function precision) */
                if (arb_is_positive(t))
                {
                    flint_printf("FAIL (%s): floor bound invalid\n", f->name);
                    goto fail;
                }
                arf_set_si_2exp_si(ulp, -3, NFLOAT_EXP(z) - fprec);
                if (arf_cmp(arb_midref(t), ulp) < 0)
                {
                    flint_printf("FAIL (%s): floor bound not tight\n", f->name);
                    goto fail;
                }
            }
            else if (flags & NFLOAT_RND_CEIL)
            {
                if (arb_is_negative(t))
                {
                    flint_printf("FAIL (%s): ceil bound invalid\n", f->name);
                    goto fail;
                }
                arf_set_si_2exp_si(ulp, 3, NFLOAT_EXP(z) - fprec);
                if (arf_cmp(arb_midref(t), ulp) > 0)
                {
                    flint_printf("FAIL (%s): ceil bound not tight\n", f->name);
                    goto fail;
                }
            }
            else
            {
                /* |z - ref| < ulp(z) at the function precision, with a small
                   tolerance */
                mag_t m, u;
                mag_init(m);
                mag_init(u);
                arb_get_mag(m, t);
                mag_set_ui_2exp_si(u, 1025, NFLOAT_EXP(z) - fprec - 10);
                if (mag_cmp(m, u) > 0)
                {
                    flint_printf("FAIL (%s): inaccurate\n", f->name);
                    mag_clear(m);
                    mag_clear(u);
                    goto fail;
                }
                mag_clear(m);
                mag_clear(u);
            }
        }
        else if (status == GR_SUCCESS)
        {
            /* special result (zero): must be a valid value */
            if (!NFLOAT_IS_ZERO(z) && !(flags & (NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN)))
            {
                flint_printf("FAIL (%s): special value\n", f->name);
                goto fail;
            }
            if (NFLOAT_IS_ZERO(z))
            {
                if (!reference(ref, f, x, (f->arity == 2) ? y : NULL, -WORD_MAX / 4, ctx) || !arb_is_zero(ref))
                {
                    /* allowed for underflow only */
                    if (arb_is_finite(ref) && !arb_is_zero(ref) &&
                        arf_cmpabs_2exp_si(arb_midref(ref), NFLOAT_MIN_EXP + 2) > 0)
                    {
                        flint_printf("FAIL (%s): zero\n", f->name);
                        goto fail;
                    }
                }
            }
        }
        else
        {
            /* failure is acceptable for overflow, underflow and domain
               errors */
            if (reference(ref, f, x, (f->arity == 2) ? y : NULL, WORD_MAX / 4, ctx) && arb_rel_accuracy_bits(ref) > 10)
            {
                if (!arb_is_zero(ref) &&
                    arf_cmpabs_2exp_si(arb_midref(ref), NFLOAT_MAX_EXP - 2) < 0 &&
                    arf_cmpabs_2exp_si(arb_midref(ref), NFLOAT_MIN_EXP + 2) > 0)
                {
                    flint_printf("FAIL (%s): unexpected status %d\n", f->name, status);
                    goto fail;
                }
            }
        }

        arb_clear(ref);
        arb_clear(t);
        arf_clear(zf);
        arf_clear(ulp);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_ctx_clear(ctx);
        continue;

fail:
        flint_printf("prec = %wd, func_prec = %wd, flags = %d, aliasing = %d\n", prec, fprec, flags, aliasing);
        flint_printf("x = %{gr}\n", x, ctx);
        flint_printf("x exp = %wd\n", NFLOAT_EXP(x));
        if (f->arity == 2)
            flint_printf("y = %{gr}\n", y, ctx);
        flint_printf("z = %{gr}\n", z, ctx);
        flint_printf("ref = "); arb_printn(ref, prec / 3.32 + 10, 0); flint_printf("\n");
        flint_printf("z - ref = "); arb_printn(t, 10, 0); flint_printf("\n");
        flint_abort();
    }

    /* the functions with two outputs at special values, with either output
       aliased to the input */
    {
        typedef int (*func2_t)(nfloat_ptr, nfloat_ptr, nfloat_srcptr, gr_ctx_t);
        static const func2_t funcs2[] = { nfloat_sin_cos, nfloat_sin_cos_pi, nfloat_sinh_cosh };
        gr_ctx_t ctx;
        nfloat_ptr x, a, b, u, w;
        int fi, v, al;

        nfloat_ctx_init(ctx, FLINT_BITS * 2, NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN);
        x = gr_heap_init(ctx);
        a = gr_heap_init(ctx);
        b = gr_heap_init(ctx);
        u = gr_heap_init(ctx);
        w = gr_heap_init(ctx);

        for (fi = 0; fi < 3; fi++)
        {
            for (v = 0; v < 4; v++)
            {
                for (al = 1; al <= 2; al++)
                {
                    int s1, s2, same;

                    if (v == 0) GR_MUST_SUCCEED(nfloat_zero(x, ctx));
                    else if (v == 1) GR_MUST_SUCCEED(nfloat_pos_inf(x, ctx));
                    else if (v == 2) GR_MUST_SUCCEED(nfloat_neg_inf(x, ctx));
                    else GR_MUST_SUCCEED(nfloat_nan(x, ctx));

                    /* (a, b) without aliasing, (u, w) with u resp. w = x */
                    s1 = funcs2[fi](a, b, x, ctx);
                    if (al == 1)
                    {
                        GR_MUST_SUCCEED(nfloat_set(u, x, ctx));
                        s2 = funcs2[fi](u, w, u, ctx);
                    }
                    else
                    {
                        GR_MUST_SUCCEED(nfloat_set(w, x, ctx));
                        s2 = funcs2[fi](u, w, w, ctx);
                    }

                    same = ((NFLOAT_IS_NAN(a) && NFLOAT_IS_NAN(u)) || gr_equal(a, u, ctx) == T_TRUE)
                        && ((NFLOAT_IS_NAN(b) && NFLOAT_IS_NAN(w)) || gr_equal(b, w, ctx) == T_TRUE);

                    if (s1 != s2 || !same)
                    {
                        flint_printf("FAIL (two outputs, aliasing): func %d, value %d, aliased output %d\n", fi, v, al);
                        flint_printf("want %{gr}, %{gr}\n", a, ctx, b, ctx);
                        flint_printf("got  %{gr}, %{gr}\n", u, ctx, w, ctx);
                        flint_abort();
                    }
                }
            }
        }

        gr_heap_clear(x, ctx);
        gr_heap_clear(a, ctx);
        gr_heap_clear(b, ctx);
        gr_heap_clear(u, ctx);
        gr_heap_clear(w, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

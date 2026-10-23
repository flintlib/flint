/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "gr.h"
#include "gr_special.h"
#include "nfloat.h"

/* Checks the complex elementary functions against acb: the result must
   lie within a few ulp (relative to the magnitude of the exact value) of
   the image of a ball around the argument with a relative radius of a
   few ulp (accounting for the conditioning of the function). */

/* (nfloat_complex_get_acb is not declared when nfloat.h comes before
   acb.h in the combined test program) */
static int
_test_get_acb(acb_t res, nfloat_complex_srcptr x, gr_ctx_t ctx)
{
    int status;
    status = nfloat_get_arf(arb_midref(acb_realref(res)), NFLOAT_COMPLEX_RE(x, ctx), ctx);
    status |= nfloat_get_arf(arb_midref(acb_imagref(res)), NFLOAT_COMPLEX_IM(x, ctx), ctx);
    mag_zero(arb_radref(acb_realref(res)));
    mag_zero(arb_radref(acb_imagref(res)));
    return status;
}

typedef struct
{
    const char * name;
    int method;
    int arity;
}
celem_func_t;

static const celem_func_t celem_funcs[] =
{
    { "exp", GR_METHOD_EXP, 1 },
    { "expm1", GR_METHOD_EXPM1, 1 },
    { "exp_pi_i", GR_METHOD_EXP_PI_I, 1 },
    { "exp2", GR_METHOD_EXP2, 1 },
    { "exp10", GR_METHOD_EXP10, 1 },
    { "log", GR_METHOD_LOG, 1 },
    { "log1p", GR_METHOD_LOG1P, 1 },
    { "log_pi_i", GR_METHOD_LOG_PI_I, 1 },
    { "log2", GR_METHOD_LOG2, 1 },
    { "log10", GR_METHOD_LOG10, 1 },
    { "sin", GR_METHOD_SIN, 1 },
    { "cos", GR_METHOD_COS, 1 },
    { "tan", GR_METHOD_TAN, 1 },
    { "cot", GR_METHOD_COT, 1 },
    { "sec", GR_METHOD_SEC, 1 },
    { "csc", GR_METHOD_CSC, 1 },
    { "sin_pi", GR_METHOD_SIN_PI, 1 },
    { "cos_pi", GR_METHOD_COS_PI, 1 },
    { "tan_pi", GR_METHOD_TAN_PI, 1 },
    { "cot_pi", GR_METHOD_COT_PI, 1 },
    { "sec_pi", GR_METHOD_SEC_PI, 1 },
    { "csc_pi", GR_METHOD_CSC_PI, 1 },
    { "sinc", GR_METHOD_SINC, 1 },
    { "sinc_pi", GR_METHOD_SINC_PI, 1 },
    { "sinh", GR_METHOD_SINH, 1 },
    { "cosh", GR_METHOD_COSH, 1 },
    { "tanh", GR_METHOD_TANH, 1 },
    { "coth", GR_METHOD_COTH, 1 },
    { "sech", GR_METHOD_SECH, 1 },
    { "csch", GR_METHOD_CSCH, 1 },
    { "asin", GR_METHOD_ASIN, 1 },
    { "acos", GR_METHOD_ACOS, 1 },
    { "atan", GR_METHOD_ATAN, 1 },
    { "acot", GR_METHOD_ACOT, 1 },
    { "asec", GR_METHOD_ASEC, 1 },
    { "acsc", GR_METHOD_ACSC, 1 },
    { "asinh", GR_METHOD_ASINH, 1 },
    { "acosh", GR_METHOD_ACOSH, 1 },
    { "atanh", GR_METHOD_ATANH, 1 },
    { "acoth", GR_METHOD_ACOTH, 1 },
    { "asech", GR_METHOD_ASECH, 1 },
    { "acsch", GR_METHOD_ACSCH, 1 },
    { "asin_pi", GR_METHOD_ASIN_PI, 1 },
    { "acos_pi", GR_METHOD_ACOS_PI, 1 },
    { "atan_pi", GR_METHOD_ATAN_PI, 1 },
    { "acot_pi", GR_METHOD_ACOT_PI, 1 },
    { "asec_pi", GR_METHOD_ASEC_PI, 1 },
    { "acsc_pi", GR_METHOD_ACSC_PI, 1 },
    { "pow", GR_METHOD_POW, 2 },
};

#define NUM_CFUNCS (sizeof(celem_funcs) / sizeof(celem_func_t))

static int
call_cfunc(const celem_func_t * f, gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    if (f->arity == 1)
        return ((gr_method_unary_op) ctx->methods[f->method])(res, x, ctx);
    else
        return ((gr_method_binary_op) ctx->methods[f->method])(res, x, y, ctx);
}

static void
randtest_real(nfloat_ptr x, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, n = NFLOAT_CTX_NLIMBS(ctx);
    slong r = n_randint(state, 10);

    if (r == 0)
    {
        nfloat_zero(x, ctx);
        return;
    }

    if (r == 1)
    {
        GR_MUST_SUCCEED(nfloat_set_si(x, (slong) n_randint(state, 7) - 3, ctx));
        if (!NFLOAT_IS_ZERO(x) && n_randint(state, 2))
            NFLOAT_EXP(x) -= 1;
        return;
    }

    for (i = 0; i < n; i++)
        NFLOAT_D(x)[i] = n_randtest(state);
    NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    NFLOAT_SGNBIT(x) = n_randint(state, 2);

    if (r == 2)
    {
        /* few bits */
        for (i = 0; i < n - 1; i++)
            NFLOAT_D(x)[i] = 0;
        NFLOAT_D(x)[n - 1] &= ~((UWORD(1) << n_randint(state, FLINT_BITS)) - 1);
        NFLOAT_D(x)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
        NFLOAT_EXP(x) = (slong) n_randint(state, 10) - 4;
    }
    else if (r == 3)
        NFLOAT_EXP(x) = (slong) n_randint(state, 400) - 200;
    else
        NFLOAT_EXP(x) = (slong) n_randint(state, 12) - 6;
}

/* acb reference for x (a ball of relative radius 2^-rad_bits); returns 0
   on failure */
static int
creference(acb_t ref, const celem_func_t * f, nfloat_complex_srcptr x, nfloat_complex_srcptr y, slong rad_bits, gr_ctx_t ctx)
{
    gr_ctx_t actx;
    acb_t ax, ay;
    mag_t r;
    slong prec = NFLOAT_CTX_PREC(ctx), wp;
    int status, ok = 0;

    acb_init(ax);
    acb_init(ay);
    mag_init(r);

    GR_MUST_SUCCEED(_test_get_acb(ax, x, ctx));
    if (rad_bits != 0)
    {
        acb_get_mag(r, ax);
        mag_mul_2exp_si(r, r, -rad_bits);
        acb_add_error_mag(ax, r);
    }
    if (y != NULL)
    {
        GR_MUST_SUCCEED(_test_get_acb(ay, y, ctx));
        if (rad_bits != 0)
        {
            acb_get_mag(r, ay);
            mag_mul_2exp_si(r, r, -rad_bits);
            acb_add_error_mag(ay, r);
        }
    }

    for (wp = 2 * prec + 64; wp < 16 * prec + 2000; wp *= 2)
    {
        gr_ctx_init_complex_acb(actx, wp);
        status = call_cfunc(f, ref, ax, ay, actx);
        gr_ctx_clear(actx);

        if (status != GR_SUCCESS || !acb_is_finite(ref))
            break;

        if (rad_bits != 0 || acb_rel_accuracy_bits(ref) > prec + 30)
        {
            ok = 1;
            break;
        }
    }

    acb_clear(ax);
    acb_clear(ay);
    mag_clear(r);
    return ok;
}

TEST_FUNCTION_START(nfloat_complex_elem, state)
{
    slong iter, count_unable = 0, count_verified = 0, count_skipped = 0;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong n, prec, fprec;
        const celem_func_t * f;
        nfloat_complex_ptr x, y, z;
        acb_t ref, ref2, t;
        mag_t m, tol;
        int status;

        if (n_randint(state, 4) == 0)
            n = 1 + n_randint(state, NFLOAT_MAX_LIMBS);
        else
            n = 1 + n_randint(state, 4);

        prec = n * FLINT_BITS;
        GR_MUST_SUCCEED(nfloat_complex_ctx_init(ctx, prec, 0));

        fprec = prec;
        if (n_randint(state, 4) == 0)
        {
            fprec = 20 + n_randint(state, prec - 19);
            GR_MUST_SUCCEED(nfloat_ctx_set_func_prec(ctx, fprec));
        }

        f = celem_funcs + n_randint(state, NUM_CFUNCS);

        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        acb_init(ref);
        acb_init(ref2);
        acb_init(t);
        mag_init(m);
        mag_init(tol);

        randtest_real(NFLOAT_COMPLEX_RE(x, ctx), state, ctx);
        randtest_real(NFLOAT_COMPLEX_IM(x, ctx), state, ctx);
        randtest_real(NFLOAT_COMPLEX_RE(y, ctx), state, ctx);
        randtest_real(NFLOAT_COMPLEX_IM(y, ctx), state, ctx);
        if (f->arity == 2 && n_randint(state, 2))
        {
            NFLOAT_EXP(NFLOAT_COMPLEX_RE(y, ctx)) = FLINT_MIN(NFLOAT_EXP(NFLOAT_COMPLEX_RE(y, ctx)), 8);
            nfloat_zero(NFLOAT_COMPLEX_IM(y, ctx), ctx);
        }

        if (n_randint(state, 2))
        {
            GR_MUST_SUCCEED(gr_set(z, x, ctx));
            status = call_cfunc(f, z, z, y, ctx);
        }
        else
            status = call_cfunc(f, z, x, y, ctx);

        if (status == GR_SUCCESS)
        {
            /* the exact value at the point (for the magnitude) and the
               image of a ball */
            /* (the acb compositions for some functions lose all accuracy
               for huge arguments) */
            if (!creference(ref, f, x, (f->arity == 2) ? y : NULL, 0, ctx))
            {
                count_skipped++;
                goto cleanup;
            }

            if (!creference(ref2, f, x, (f->arity == 2) ? y : NULL, fprec - 10, ctx))
            {
                /* (the ball contains a singularity) */
                count_skipped++;
                goto cleanup;
            }

            /* z must lie in ref2 + 2^-(fprec - 10) |ref| */
            GR_MUST_SUCCEED(_test_get_acb(t, z, ctx));
            acb_get_mag(tol, ref);
            mag_mul_2exp_si(tol, tol, -(fprec - 10));
            acb_add_error_mag(ref2, tol);

            if (!acb_contains(ref2, t))
            {
                flint_printf("FAIL (%s): inaccurate\n", f->name);
                goto fail;
            }

            count_verified++;
        }
        else
        {
            count_unable++;

            /* failure is acceptable near singularities and for huge or
               tiny results */
            if (creference(ref, f, x, (f->arity == 2) ? y : NULL, 0, ctx) &&
                creference(ref2, f, x, (f->arity == 2) ? y : NULL, fprec - 8, ctx))
            {
                arb_srcptr re = acb_realref(ref), im = acb_imagref(ref);

                /* (each nonzero component must be in the nfloat range) */
                acb_get_mag(m, ref);
                if (mag_cmp_2exp_si(m, 100000) < 0 && mag_cmp_2exp_si(m, -100000) > 0 &&
                    (arb_is_zero(re) || arf_cmpabs_2exp_si(arb_midref(re), NFLOAT_MIN_EXP + 2) > 0) &&
                    (arb_is_zero(im) || arf_cmpabs_2exp_si(arb_midref(im), NFLOAT_MIN_EXP + 2) > 0) &&
                    (acb_is_zero(ref) || acb_rel_accuracy_bits(ref2) > 10))
                {
                    flint_printf("FAIL (%s): unexpected status %d\n", f->name, status);
                    goto fail;
                }
            }
        }

cleanup:
        acb_clear(ref);
        acb_clear(ref2);
        acb_clear(t);
        mag_clear(m);
        mag_clear(tol);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_ctx_clear(ctx);
        continue;

fail:
        flint_printf("prec = %wd, func_prec = %wd\n", prec, fprec);
        flint_printf("x = %{gr}\n", x, ctx);
        if (f->arity == 2)
            flint_printf("y = %{gr}\n", y, ctx);
        flint_printf("z = %{gr}\n", z, ctx);
        flint_printf("ref = "); acb_printn(ref, 30, 0); flint_printf("\n");
        flint_printf("ref2 = "); acb_printn(ref2, 30, 0); flint_printf("\n");
        flint_abort();
    }

    if (getenv("NFLOAT_TEST_VERBOSE"))
        flint_printf("verified %wd, skipped %wd, unable %wd\n", count_verified, count_skipped, count_unable);

    TEST_FUNCTION_END(state);
}

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Timings (ns per call) of the nfloat elementary functions against
   MPFR and arb at the same precision.

   Usage: p-elem [func [flags [func_prec_deficit]]]

   func is a function name (default: all), flags the nfloat context
   flags (8 = NFLOAT_RND_FLOOR), and func_prec_deficit a number of bits
   subtracted from the function precision. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpfr.h>
#include "arb.h"
#include "gr.h"
#include "gr_special.h"
#include "nfloat.h"
#include "profiler.h"

#define N 128

typedef int (*mpfr_func1)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
typedef int (*mpfr_func2)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);

typedef struct
{
    const char * name;
    int method;
    int arity;
    void (* mpfr_func)(void);   /* an mpfr_func1 or mpfr_func2 */
    double lo, hi;          /* argument range: |x| in 2^[lo, hi] */
    int positive;
}
func_t;

static int _mpfr_sin(mpfr_ptr r, mpfr_srcptr x, mpfr_rnd_t rnd) { return mpfr_sin(r, x, rnd); }

static const func_t funcs[] =
{
    { "exp", GR_METHOD_EXP, 1, (void (*)(void)) mpfr_exp, -4, 4, 0 },
    { "expm1", GR_METHOD_EXPM1, 1, (void (*)(void)) mpfr_expm1, -4, 4, 0 },
    { "exp2", GR_METHOD_EXP2, 1, (void (*)(void)) mpfr_exp2, -4, 4, 0 },
    { "log", GR_METHOD_LOG, 1, (void (*)(void)) mpfr_log, -4, 4, 1 },
    { "log1p", GR_METHOD_LOG1P, 1, (void (*)(void)) mpfr_log1p, -4, 4, 1 },
    { "log2", GR_METHOD_LOG2, 1, (void (*)(void)) mpfr_log2, -4, 4, 1 },
    { "sin", GR_METHOD_SIN, 1, (void (*)(void)) _mpfr_sin, -4, 4, 0 },
    { "cos", GR_METHOD_COS, 1, (void (*)(void)) mpfr_cos, -4, 4, 0 },
    { "tan", GR_METHOD_TAN, 1, (void (*)(void)) mpfr_tan, -4, 4, 0 },
    { "atan", GR_METHOD_ATAN, 1, (void (*)(void)) mpfr_atan, -4, 4, 0 },
    { "atan2", GR_METHOD_ATAN2, 2, (void (*)(void)) mpfr_atan2, -4, 4, 0 },
    { "asin", GR_METHOD_ASIN, 1, (void (*)(void)) mpfr_asin, -4, -0.01, 0 },
    { "acos", GR_METHOD_ACOS, 1, (void (*)(void)) mpfr_acos, -4, -0.01, 0 },
    { "sinh", GR_METHOD_SINH, 1, (void (*)(void)) mpfr_sinh, -4, 4, 0 },
    { "cosh", GR_METHOD_COSH, 1, (void (*)(void)) mpfr_cosh, -4, 4, 0 },
    { "tanh", GR_METHOD_TANH, 1, (void (*)(void)) mpfr_tanh, -4, 4, 0 },
    { "asinh", GR_METHOD_ASINH, 1, (void (*)(void)) mpfr_asinh, -4, 4, 0 },
    { "acosh", GR_METHOD_ACOSH, 1, (void (*)(void)) mpfr_acosh, 0.01, 4, 1 },
    { "atanh", GR_METHOD_ATANH, 1, (void (*)(void)) mpfr_atanh, -4, -0.01, 0 },
    { "sin_pi", GR_METHOD_SIN_PI, 1, (void (*)(void)) mpfr_sinpi, -4, 4, 0 },
    { "pow", GR_METHOD_POW, 2, (void (*)(void)) mpfr_pow, -2, 2, 1 },
};

#define NUM_FUNCS (sizeof(funcs) / sizeof(func_t))

static double
now(void)
{
    timeit_t t;
    timeit_start_us(t);
    return -1e-6 * (double) t->wall;
}

/* best of 5 runs of body over the N arguments, in ns per call */
#define TIME(result, ...) \
    do { \
        double __best = 1e300; \
        int __k; \
        for (__k = 0; __k < 5; __k++) \
        { \
            slong __reps = 1, __r; \
            double __t; \
            for (;;) \
            { \
                __t = now(); \
                for (__r = 0; __r < __reps; __r++) { __VA_ARGS__; } \
                __t = now() - __t; \
                if (__t > 0.01) break; \
                __reps *= 2; \
            } \
            __t = __t / __reps; \
            if (__t < __best) __best = __t; \
        } \
        (result) = __best * 1e9 / N; \
    } while (0)

int main(int argc, char * argv[])
{
    slong nlist[] = { 1, 2, 3, 4, 6, 8, 16, 32, 66, 0 };
    const char * only = (argc > 1) ? argv[1] : NULL;
    int flags = (argc > 2) ? atoi(argv[2]) : 0;
    slong deficit = (argc > 3) ? atol(argv[3]) : 0;
    mpfr_rnd_t rnd;
    slong fi, i, j;
    flint_rand_t state;

    if (flags & NFLOAT_RND_FLOOR)
        rnd = MPFR_RNDD;
    else if (flags & NFLOAT_RND_CEIL)
        rnd = MPFR_RNDU;
    else
        rnd = MPFR_RNDZ;

    flint_rand_init(state);

    flint_printf("ns per call; nfloat (flags %d, func_prec deficit %wd), speedup vs mpfr, arb, old arf path\n", flags, deficit);

    for (fi = 0; fi < (slong) NUM_FUNCS; fi++)
    {
        const func_t * f = funcs + fi;

        if (only != NULL && strcmp(only, "all") && strcmp(only, f->name))
            continue;

        flint_printf("%-8s", f->name);
        for (i = 0; nlist[i] != 0; i++)
            flint_printf(" | %16wd", nlist[i] * FLINT_BITS);
        flint_printf("\n");

        flint_printf("%-8s", "");
        for (i = 0; nlist[i] != 0; i++)
        {
            slong n = nlist[i], prec = n * FLINT_BITS, sz;
            gr_ctx_t ctx, actx, rctx;
            gr_ptr X, Y, Z;
            mpfr_ptr fx, fy, fz;
            arb_ptr ax, ay, az;
            gr_method_unary_op op1;
            gr_method_binary_op op2;
            int status = GR_SUCCESS;
            double tn, tm, ta, to;

            nfloat_ctx_init(ctx, prec, flags);
            if (deficit)
                GR_MUST_SUCCEED(nfloat_ctx_set_func_prec(ctx, prec - deficit));
            gr_ctx_init_real_arb(actx, prec);
            gr_ctx_init_real_float_arf(rctx, prec);
            sz = ctx->sizeof_elem;

            X = gr_heap_init_vec(N, ctx);
            Y = gr_heap_init_vec(N, ctx);
            Z = gr_heap_init_vec(N, ctx);
            fx = flint_malloc(N * sizeof(__mpfr_struct));
            fy = flint_malloc(N * sizeof(__mpfr_struct));
            fz = flint_malloc(N * sizeof(__mpfr_struct));
            ax = _arb_vec_init(N);
            ay = _arb_vec_init(N);
            az = _arb_vec_init(N);

            for (j = 0; j < N; j++)
            {
                arf_t t;
                double d;

                d = f->lo + (f->hi - f->lo) * (n_randlimb(state) % 1000) / 1000.0;
                d = pow(2.0, d);
                if (!f->positive && n_randint(state, 2))
                    d = -d;

                arf_init(t);
                arf_set_d(t, d);
                /* fill all mantissa bits */
                {
                    arf_t u;
                    arf_init(u);
                    arf_urandom(u, state, prec, ARF_RND_DOWN);
                    arf_mul_2exp_si(u, u, -60);
                    arf_mul(u, u, t, ARF_PREC_EXACT, ARF_RND_DOWN);
                    arf_add(t, t, u, prec, ARF_RND_DOWN);
                    arf_clear(u);
                }
                GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(X, j, sz), t, ctx));
                GR_MUST_SUCCEED(nfloat_get_arf(arb_midref(ax + j), GR_ENTRY(X, j, sz), ctx));
                mpfr_init2(fx + j, prec);
                mpfr_init2(fy + j, prec);
                mpfr_init2(fz + j, prec);
                arf_get_mpfr(fx + j, arb_midref(ax + j), MPFR_RNDN);

                d = pow(2.0, f->lo + (f->hi - f->lo) * (n_randlimb(state) % 1000) / 1000.0);
                arf_set_d(t, d);
                GR_MUST_SUCCEED(nfloat_set_arf(GR_ENTRY(Y, j, sz), t, ctx));
                GR_MUST_SUCCEED(nfloat_get_arf(arb_midref(ay + j), GR_ENTRY(Y, j, sz), ctx));
                arf_get_mpfr(fy + j, arb_midref(ay + j), MPFR_RNDN);
                arf_clear(t);
            }

            if (f->arity == 1)
            {
                op1 = (gr_method_unary_op) ctx->methods[f->method];
                TIME(tn, for (j = 0; j < N; j++) status |= op1(GR_ENTRY(Z, j, sz), GR_ENTRY(X, j, sz), ctx));
                TIME(tm, for (j = 0; j < N; j++) ((mpfr_func1) f->mpfr_func)(fz + j, fx + j, rnd));
                op1 = (gr_method_unary_op) actx->methods[f->method];
                TIME(ta, for (j = 0; j < N; j++) status |= op1(az + j, ax + j, actx));
                op1 = (gr_method_unary_op) rctx->methods[f->method];
                TIME(to, for (j = 0; j < N; j++) { arf_t t; arf_init(t); nfloat_get_arf(t, GR_ENTRY(X, j, sz), ctx); op1(t, t, rctx); nfloat_set_arf(GR_ENTRY(Z, j, sz), t, ctx); arf_clear(t); });
            }
            else
            {
                op2 = (gr_method_binary_op) ctx->methods[f->method];
                TIME(tn, for (j = 0; j < N; j++) status |= op2(GR_ENTRY(Z, j, sz), GR_ENTRY(X, j, sz), GR_ENTRY(Y, j, sz), ctx));
                TIME(tm, for (j = 0; j < N; j++) ((mpfr_func2) f->mpfr_func)(fz + j, fx + j, fy + j, rnd));
                op2 = (gr_method_binary_op) actx->methods[f->method];
                TIME(ta, for (j = 0; j < N; j++) status |= op2(az + j, ax + j, ay + j, actx));
                op2 = (gr_method_binary_op) rctx->methods[f->method];
                TIME(to, for (j = 0; j < N; j++) { arf_t t, u; arf_init(t); arf_init(u); nfloat_get_arf(t, GR_ENTRY(X, j, sz), ctx); nfloat_get_arf(u, GR_ENTRY(Y, j, sz), ctx); op2(t, t, u, rctx); nfloat_set_arf(GR_ENTRY(Z, j, sz), t, ctx); arf_clear(t); arf_clear(u); });
            }

            flint_printf(" | %5.0f %3.0f %3.1f %3.1f", tn, tm / tn, ta / tn, to / tn);
            fflush(stdout);
            (void) status;

            for (j = 0; j < N; j++)
            {
                mpfr_clear(fx + j);
                mpfr_clear(fy + j);
                mpfr_clear(fz + j);
            }
            flint_free(fx);
            flint_free(fy);
            flint_free(fz);
            _arb_vec_clear(ax, N);
            _arb_vec_clear(ay, N);
            _arb_vec_clear(az, N);
            gr_heap_clear_vec(X, N, ctx);
            gr_heap_clear_vec(Y, N, ctx);
            gr_heap_clear_vec(Z, N, ctx);
            gr_ctx_clear(ctx);
            gr_ctx_clear(actx);
            gr_ctx_clear(rctx);
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

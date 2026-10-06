/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* ns per operation for dN, dNb (with an s suffix: the strong,
   canonicalising contexts; f: the fast contexts that bound instead of
   track the rounding errors of products), nfloat and arb at the same
   precision, on random inputs of moderate magnitude, in two tables:
   the arithmetic (scalar operations timed as independent operations
   over an array, i.e. throughput, and the vector methods add, mul and
   dot per element), and the elementary functions, each as a scalar
   method and as a vector function per element (dN, dNb and arb). */

#include <stdlib.h>
#include "profiler.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "arb.h"
#include "nfloat.h"
#include "dfloat.h"

#define LEN 1024

static void
random_vec(gr_ptr x, slong len, flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        /* values near 1 with full precision */
        GR_MUST_SUCCEED(gr_set_ui(GR_ENTRY(x, i, ctx->sizeof_elem), 1 + n_randint(state, 1000), ctx));
        GR_MUST_SUCCEED(gr_sqrt(GR_ENTRY(x, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx));
    }
}

#define TIME_BEST(best, body) \
    do { \
        int _pass; slong _reps; timeit_t _tm; \
        (best) = 1e300; \
        for (_pass = 0; _pass < 5; _pass++) \
        { \
            _reps = 1; \
            while (1) \
            { \
                slong _r; \
                timeit_start(_tm); \
                for (_r = 0; _r < _reps; _r++) { body; } \
                timeit_stop(_tm); \
                if (_tm->wall >= 20) break; \
                _reps *= 4; \
            } \
            (best) = FLINT_MIN((best), 1e6 * _tm->wall / (_reps * (double) LEN)); \
        } \
    } while (0)

static double
time_op(gr_ctx_t ctx, gr_method method, int unary)
{
    gr_ptr x, y, z;
    slong i;
    double best;
    flint_rand_t state;

    flint_rand_init(state);
    x = gr_heap_init_vec(LEN, ctx);
    y = gr_heap_init_vec(LEN, ctx);
    z = gr_heap_init_vec(LEN, ctx);
    random_vec(x, LEN, state, ctx);
    random_vec(y, LEN, state, ctx);

    TIME_BEST(best,
        for (i = 0; i < LEN; i++)
        {
            if (unary)
                GR_IGNORE(((gr_method_unary_op *) ctx->methods)[method](GR_ENTRY(z, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx));
            else
                GR_IGNORE(((gr_method_binary_op *) ctx->methods)[method](GR_ENTRY(z, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), GR_ENTRY(y, i, ctx->sizeof_elem), ctx));
        });

    gr_heap_clear_vec(x, LEN, ctx);
    gr_heap_clear_vec(y, LEN, ctx);
    gr_heap_clear_vec(z, LEN, ctx);
    flint_rand_clear(state);
    return best;
}

/* which: 0 = vec_add, 1 = vec_mul, 2 = vec_dot */
static double
time_vec(gr_ctx_t ctx, int which)
{
    gr_ptr x, y, z;
    double best;
    flint_rand_t state;

    flint_rand_init(state);
    x = gr_heap_init_vec(LEN, ctx);
    y = gr_heap_init_vec(LEN, ctx);
    z = gr_heap_init_vec(LEN, ctx);
    random_vec(x, LEN, state, ctx);
    random_vec(y, LEN, state, ctx);

    if (which == 0)
        TIME_BEST(best, GR_IGNORE(_gr_vec_add(z, x, y, LEN, ctx)));
    else if (which == 1)
        TIME_BEST(best, GR_IGNORE(_gr_vec_mul(z, x, y, LEN, ctx)));
    else
        TIME_BEST(best, GR_IGNORE(_gr_vec_dot(z, NULL, 0, x, y, LEN, ctx)));

    gr_heap_clear_vec(x, LEN, ctx);
    gr_heap_clear_vec(y, LEN, ctx);
    gr_heap_clear_vec(z, LEN, ctx);
    flint_rand_clear(state);
    return best;
}

/* the vector elementary functions: the typed ones for dfloat (all of
   them have one, atan included), else the generic ring (exp, sin, log)
   or a loop of the scalar method (atan) */
#define VEC_TYPED(n, B, f) \
    do { \
        if (n == 1) GR_IGNORE((_d1##B##_vec_##f((void *) z, (const void *) x, LEN), GR_SUCCESS)); \
        else if (n == 2) GR_IGNORE((_d2##B##_vec_##f((void *) z, (const void *) x, LEN), GR_SUCCESS)); \
        else if (n == 3) GR_IGNORE((_d3##B##_vec_##f((void *) z, (const void *) x, LEN), GR_SUCCESS)); \
        else GR_IGNORE((_d4##B##_vec_##f((void *) z, (const void *) x, LEN), GR_SUCCESS)); \
    } while (0)

/* which: 0 = exp, 1 = sin, 2 = log, 3 = atan */
static double
time_vec_elem(gr_ctx_t ctx, int which, int dfloat)
{
    gr_ptr x, z;
    double best;
    flint_rand_t state;
    slong i;

    flint_rand_init(state);
    x = gr_heap_init_vec(LEN, ctx);
    z = gr_heap_init_vec(LEN, ctx);
    random_vec(x, LEN, state, ctx);

    if (dfloat)
    {
        int n = DFLOAT_CTX_N(ctx), ball = DFLOAT_CTX_BALL(ctx);
        if (which == 0) { if (ball) TIME_BEST(best, VEC_TYPED(n, b, exp)); else TIME_BEST(best, VEC_TYPED(n, , exp)); }
        else if (which == 1) { if (ball) TIME_BEST(best, VEC_TYPED(n, b, sin)); else TIME_BEST(best, VEC_TYPED(n, , sin)); }
        else if (which == 2) { if (ball) TIME_BEST(best, VEC_TYPED(n, b, log)); else TIME_BEST(best, VEC_TYPED(n, , log)); }
        else { if (ball) TIME_BEST(best, VEC_TYPED(n, b, atan)); else TIME_BEST(best, VEC_TYPED(n, , atan)); }
    }
    else if (which == 0)
        TIME_BEST(best, GR_IGNORE(_gr_vec_exp(z, x, LEN, ctx)));
    else if (which == 1)
        TIME_BEST(best, GR_IGNORE(_gr_vec_sin(z, x, LEN, ctx)));
    else if (which == 2)
        TIME_BEST(best, GR_IGNORE(_gr_vec_log(z, x, LEN, ctx)));
    else
        TIME_BEST(best,
            for (i = 0; i < LEN; i++)
                GR_IGNORE(gr_atan(GR_ENTRY(z, i, ctx->sizeof_elem), GR_ENTRY(x, i, ctx->sizeof_elem), ctx)));

    gr_heap_clear_vec(x, LEN, ctx);
    gr_heap_clear_vec(z, LEN, ctx);
    flint_rand_clear(state);
    return best;
}

/* the arithmetic: scalar operations, then the vector ones */
static void
row_arith(const char * name, gr_ctx_t ctx)
{
    const gr_method methods[] = { GR_METHOD_ADD, GR_METHOD_MUL, GR_METHOD_SQR, GR_METHOD_DIV, GR_METHOD_SQRT };
    const int unary[] = { 0, 0, 1, 0, 1 };
    int k;

    flint_printf("%-10s", name);
    for (k = 0; k < 5; k++)
        flint_printf("%8.1f", time_op(ctx, methods[k], unary[k]));
    flint_printf("  ");
    for (k = 0; k < 3; k++)
        flint_printf("%8.2f", time_vec(ctx, k));
    flint_printf("\n");
    fflush(stdout);
}

/* the elementary functions, each scalar and vector */
static void
row_elem(const char * name, gr_ctx_t ctx, int dfloat)
{
    const gr_method methods[] = { GR_METHOD_EXP, GR_METHOD_SIN, GR_METHOD_LOG, GR_METHOD_ATAN };
    int k;

    flint_printf("%-10s", name);
    for (k = 0; k < 4; k++)
        flint_printf("%8.1f%8.1f", time_op(ctx, methods[k], 1), time_vec_elem(ctx, k, dfloat));
    flint_printf("\n");
    fflush(stdout);
}

int main(void)
{
    gr_ctx_t ctx;
    int n, ball;
    char name[16];

    flint_printf("Arithmetic (ns per operation / element)\n");
    flint_printf("%-10s%8s%8s%8s%8s%8s  %8s%8s%8s\n",
        "ring", "add", "mul", "sqr", "div", "sqrt", "vec_add", "vec_mul", "vec_dot");

    for (n = 1; n <= 4; n++)
    {
        for (ball = 0; ball <= 1; ball++)
        {
            int strong, fast;
            for (fast = 0; fast <= ball; fast++)
                for (strong = 0; strong <= 1; strong++)
                {
                    GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, (ball ? DFLOAT_BALL : 0) | (strong ? DFLOAT_STRONG : 0) | (fast ? DFLOAT_FAST : 0)));
                    flint_sprintf(name, "d%d%s%s%s", n, ball ? "b" : "", fast ? "f" : "", strong ? "s" : "");
                    row_arith(name, ctx);
                    gr_ctx_clear(ctx);
                }
        }
        GR_MUST_SUCCEED(nfloat_ctx_init(ctx, 64 * n, 0));
        flint_sprintf(name, "nfloat%d", 64 * n);
        row_arith(name, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_init_real_arb(ctx, 53 * n);
        flint_sprintf(name, "arb%d", 53 * n);
        row_arith(name, ctx);
        gr_ctx_clear(ctx);
    }

    flint_printf("\nElementary functions, scalar and vector (ns per element)\n");
    flint_printf("%-10s%8s%8s%8s%8s%8s%8s%8s%8s\n",
        "ring", "exp", "vec", "sin", "vec", "log", "vec", "atan", "vec");

    for (n = 1; n <= 4; n++)
    {
        for (ball = 0; ball <= 1; ball++)
        {
            GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, ball ? DFLOAT_BALL : 0));
            flint_sprintf(name, "d%d%s", n, ball ? "b" : "");
            row_elem(name, ctx, 1);
            gr_ctx_clear(ctx);
        }
        gr_ctx_init_real_arb(ctx, 53 * n);
        flint_sprintf(name, "arb%d", 53 * n);
        row_elem(name, ctx, 0);
        gr_ctx_clear(ctx);
    }
    return 0;
}

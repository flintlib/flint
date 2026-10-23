/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "decimal.h"
#include "arf.h"
#include "arb.h"
#include "gr.h"

#define N 64

static double
_time_decfloat(int op, slong prec, flint_rand_t state, int ball)
{
    gr_ctx_t ctx;
    gr_ptr x, y, z;
    slong i, reps;
    double t;
    timeit_t timer;

    if (ball)
        gr_ctx_init_decball(ctx, prec, 0);
    else
        gr_ctx_init_decfloat(ctx, prec, 0);

    x = gr_heap_init_vec(N, ctx);
    y = gr_heap_init_vec(N, ctx);
    z = gr_heap_init_vec(N, ctx);

    for (i = 0; i < N; i++)
    {
        /* full-precision random values with nearby exponents */
        decfloat_ptr xi = ball ? DECBALL_MIDREF((decball_ptr) GR_ENTRY(x, i, ctx->sizeof_elem)) : GR_ENTRY(x, i, ctx->sizeof_elem);
        decfloat_ptr yi = ball ? DECBALL_MIDREF((decball_ptr) GR_ENTRY(y, i, ctx->sizeof_elem)) : GR_ENTRY(y, i, ctx->sizeof_elem);
        fmpz_t m, e;
        fmpz_init(m);
        fmpz_init(e);
        fmpz_randbits(m, state, (slong) (prec * 3.3219) + 1);
        fmpz_abs(m, m);
        fmpz_set_si(e, -prec + n_randint(state, 3));
        GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(xi, m, e, prec, DECIMAL_RND_NEAR, ctx));
        fmpz_randbits(m, state, (slong) (prec * 3.3219) + 1);
        fmpz_abs(m, m);
        fmpz_set_si(e, -prec + n_randint(state, 3));
        GR_MUST_SUCCEED(decfloat_set_round_fmpz_10exp_fmpz(yi, m, e, prec, DECIMAL_RND_NEAR, ctx));
        fmpz_clear(m);
        fmpz_clear(e);
        if (ball)
        {
            _decmag_zero(DECBALL_RADREF((decball_ptr) GR_ENTRY(x, i, ctx->sizeof_elem)), ctx);
            _decmag_zero(DECBALL_RADREF((decball_ptr) GR_ENTRY(y, i, ctx->sizeof_elem)), ctx);
        }
    }

    reps = 1;
    for (;;)
    {
        slong r;
        timeit_start(timer);
        for (r = 0; r < reps; r++)
        {
            for (i = 0; i < N; i++)
            {
                gr_ptr xi = GR_ENTRY(x, i, ctx->sizeof_elem);
                gr_ptr yi = GR_ENTRY(y, i, ctx->sizeof_elem);
                gr_ptr zi = GR_ENTRY(z, i, ctx->sizeof_elem);
                switch (op)
                {
                    case 0: GR_IGNORE(gr_add(zi, xi, yi, ctx)); break;
                    case 1: GR_IGNORE(gr_mul(zi, xi, yi, ctx)); break;
                    case 2: GR_IGNORE(gr_div(zi, xi, yi, ctx)); break;
                    case 3: GR_IGNORE(gr_sqrt(zi, xi, ctx)); break;
                }
            }
        }
        timeit_stop(timer);
        if (timer->cpu >= 100)
            break;
        reps *= 2;
    }

    t = (double) timer->cpu / (reps * N) * 1e6;   /* ns */

    gr_heap_clear_vec(x, N, ctx);
    gr_heap_clear_vec(y, N, ctx);
    gr_heap_clear_vec(z, N, ctx);
    gr_ctx_clear(ctx);
    return t;
}

static double
_time_arf(int op, slong prec_bits, flint_rand_t state, int ball)
{
    gr_ctx_t ctx;
    gr_ptr x, y, z;
    slong i, reps;
    double t;
    timeit_t timer;

    if (ball)
        gr_ctx_init_real_arb(ctx, prec_bits);
    else
        gr_ctx_init_real_float_arf(ctx, prec_bits);

    x = gr_heap_init_vec(N, ctx);
    y = gr_heap_init_vec(N, ctx);
    z = gr_heap_init_vec(N, ctx);

    for (i = 0; i < N; i++)
    {
        arf_ptr xi = ball ? arb_midref((arb_ptr) GR_ENTRY(x, i, ctx->sizeof_elem)) : GR_ENTRY(x, i, ctx->sizeof_elem);
        arf_ptr yi = ball ? arb_midref((arb_ptr) GR_ENTRY(y, i, ctx->sizeof_elem)) : GR_ENTRY(y, i, ctx->sizeof_elem);
        arf_randtest(xi, state, prec_bits, 1);
        arf_randtest(yi, state, prec_bits, 1);
        arf_abs(xi, xi);
        arf_abs(yi, yi);
        arf_add_ui(xi, xi, 1, prec_bits, ARF_RND_DOWN);
        arf_add_ui(yi, yi, 1, prec_bits, ARF_RND_DOWN);
    }

    reps = 1;
    for (;;)
    {
        slong r;
        timeit_start(timer);
        for (r = 0; r < reps; r++)
        {
            for (i = 0; i < N; i++)
            {
                gr_ptr xi = GR_ENTRY(x, i, ctx->sizeof_elem);
                gr_ptr yi = GR_ENTRY(y, i, ctx->sizeof_elem);
                gr_ptr zi = GR_ENTRY(z, i, ctx->sizeof_elem);
                switch (op)
                {
                    case 0: GR_IGNORE(gr_add(zi, xi, yi, ctx)); break;
                    case 1: GR_IGNORE(gr_mul(zi, xi, yi, ctx)); break;
                    case 2: GR_IGNORE(gr_div(zi, xi, yi, ctx)); break;
                    case 3: GR_IGNORE(gr_sqrt(zi, xi, ctx)); break;
                }
            }
        }
        timeit_stop(timer);
        if (timer->cpu >= 100)
            break;
        reps *= 2;
    }

    t = (double) timer->cpu / (reps * N) * 1e6;

    gr_heap_clear_vec(x, N, ctx);
    gr_heap_clear_vec(y, N, ctx);
    gr_heap_clear_vec(z, N, ctx);
    gr_ctx_clear(ctx);
    return t;
}

int main(int argc, char * argv[])
{
    flint_rand_t state;
    slong precs[] = { 10, 19, 38, 100, 300, 1000, 3000, 10000, 30000, 100000, 0 };
    const char * ops[] = { "add", "mul", "div", "sqrt" };
    slong i, op;
    int ball;

    flint_rand_init(state);

    for (ball = 0; ball <= 1; ball++)
    {
        flint_printf("\n%s (ns per operation; ratio decimal/binary in parentheses)\n", ball ? "decball vs arb" : "decfloat vs arf");
        flint_printf("%10s", "digits");
        for (op = 0; op < 4; op++)
            flint_printf("  %22s", ops[op]);
        flint_printf("\n");

        for (i = 0; precs[i] != 0; i++)
        {
            slong prec = precs[i];
            slong bits = (slong) (prec * 3.3219280948873623479);

            flint_printf("%10wd", prec);
            for (op = 0; op < 4; op++)
            {
                double t1, t2;
                if (prec >= 30000 && op >= 2 && argc < 2)
                {
                    flint_printf("  %22s", "-");
                    continue;
                }
                t1 = _time_decfloat(op, prec, state, ball);
                t2 = _time_arf(op, bits, state, ball);
                flint_printf("  %12.0f (%6.2fx)", t1, t1 / t2);
                fflush(stdout);
            }
            flint_printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}

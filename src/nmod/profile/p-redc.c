/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "double_extras.h"
#include "profiler.h"

ulong pbits_arr[] = { FLINT_BITS / 2 - 2, FLINT_BITS - 8, FLINT_BITS, 0 };
ulong ebits_arr[] = { 2, 4, 8, 12, 16, 32, 48, 64, 0 };

FLINT_STATIC_NOINLINE void _loop_redc_set_nmod(nn_ptr z, nn_srcptr x, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_set_nmod(x[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_redc_get_nmod(nn_ptr z, nn_srcptr x, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_get_nmod(x[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_mul(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_t mod)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_mul(x[i], y[i], mod);
}

FLINT_STATIC_NOINLINE void _loop_nmod_mul_fullword(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_t mod)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = _nmod_mul_fullword(x[i], y[i], mod);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_mul(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_mul(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_fast_mul(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_fast_mul(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_add(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_t mod)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_add(x[i], y[i], mod);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_add(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_add(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_fast_add(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_fast_add(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_fast_mul_two(nn_ptr z, nn_srcptr x, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_fast_mul_two(x[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_redc_half_set_nmod(nn_ptr z, nn_srcptr x, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_set_nmod(x[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_redc_half_get_nmod(nn_ptr z, nn_srcptr x, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_get_nmod(x[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_half_mul(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_mul(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_half_fast_mul(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_fast_mul(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_half_add(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_add(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE void _loop_nmod_redc_half_fast_add(nn_ptr z, nn_srcptr x, nn_srcptr y, slong N, nmod_redc_ctx_t ctx)
{
    slong i;
    for (i = 0; i < N; i++)
        z[i] = nmod_redc_half_fast_add(x[i], y[i], ctx);
}

FLINT_STATIC_NOINLINE ulong
_nmod_redc_half_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;

    while ((exp & 1) == 0)
    {
        a = nmod_redc_half_mul(a, a, ctx);
        exp >>= 1;
    }

    x = a;
    while (exp >>= 1)
    {
        a = nmod_redc_half_mul(a, a, ctx);
        if (exp & 1)
            x = nmod_redc_half_mul(x, a, ctx);
    }

    return x;
}

FLINT_STATIC_NOINLINE ulong
_nmod_redc_half_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
{
    ulong x;

    while ((exp & 1) == 0)
    {
        a = nmod_redc_half_fast_mul(a, a, ctx);
        exp >>= 1;
    }

    x = a;
    while (exp >>= 1)
    {
        a = nmod_redc_half_fast_mul(a, a, ctx);
        if (exp & 1)
            x = nmod_redc_half_fast_mul(x, a, ctx);
    }

    return x;
}

int main()
{
    flint_rand_t state;
    slong N;
    ulong n;
    nmod_t mod;
    nmod_redc_ctx_t ctx, ctx_half;
    nn_ptr x, y, z, x_red, y_red, x_red_fast, y_red_fast, exps;
    nn_ptr x_red_half, y_red_half, x_red_half_fast, y_red_half_fast;
    ulong pbits, ebits;
    ulong s;
    slong i;
    slong ipbits, iebits;
    double t1, t2, t3, t4, t5, t6, FLINT_SET_BUT_UNUSED(tcpu);

    flint_rand_init(state);

    flint_printf("All times in nanoseconds\n");

    N = 1000;

    x = _nmod_vec_init(N);
    y = _nmod_vec_init(N);
    z = _nmod_vec_init(N);
    x_red = _nmod_vec_init(N);
    y_red = _nmod_vec_init(N);
    x_red_half = _nmod_vec_init(N);
    y_red_half = _nmod_vec_init(N);
    x_red_fast = _nmod_vec_init(N);
    y_red_fast = _nmod_vec_init(N);
    x_red_half_fast = _nmod_vec_init(N);
    y_red_half_fast = _nmod_vec_init(N);
    exps = _nmod_vec_init(N);

#define PRINTF(f) flint_printf("%30s    ", f)
#define PRINTF2(f, t) flint_printf("%30s    %g\n", f, t / N / 1e-9)

    s = 0;

    for (ipbits = 0; (pbits = pbits_arr[ipbits]) != 0; ipbits++)
    {
        n = n_nextprime(UWORD(1) << (pbits - 1), 1);
        nmod_init(&mod, n);
        nmod_redc_ctx_init_nmod(ctx, mod);

        if (pbits <= FLINT_BITS / 2 - 2)
            nmod_redc_ctx_init_nmod(ctx_half, mod);

        flint_printf("\nn = %wu\n\n", n);

        for (i = 0; i < N; i++)
        {
            x[i] = n_randtest(state) % n;  
            y[i] = n_randtest(state) % n;
            x_red[i] = nmod_redc_set_nmod(x[i], ctx);
            y_red[i] = nmod_redc_set_nmod(y[i], ctx);
            x_red_fast[i] = x_red[i] + n_randint(state, 2) * n;
            y_red_fast[i] = y_red[i] + n_randint(state, 2) * n;

            if (pbits <= FLINT_BITS / 2 - 2)
            {
                x_red_half[i] = nmod_redc_half_set_nmod(x[i], ctx);
                y_red_half[i] = nmod_redc_half_set_nmod(y[i], ctx);
                x_red_half_fast[i] = x_red_half[i] + n_randint(state, 2) * n;
                y_red_half_fast[i] = y_red_half[i] + n_randint(state, 2) * n;
            }
        }

        TIMEIT_START;
        for (i = 0; i < N; i++)
            nmod_redc_ctx_init_nmod(ctx, mod);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redc_ctx_init_nmod", t1);

        TIMEIT_START;
        _loop_redc_set_nmod(z, x, N, ctx);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redc_set_nmod", t1);

        TIMEIT_START;
        _loop_redc_get_nmod(z, x, N, ctx);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redc_get_nmod", t1);

        if (pbits <= FLINT_BITS / 2 - 2)
        {
            TIMEIT_START;
            for (i = 0; i < N; i++)
                nmod_redc_half_ctx_init_nmod(ctx_half, mod);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_ctx_init_nmod", t1);

            TIMEIT_START;
            _loop_redc_half_set_nmod(z, x, N, ctx_half);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_set_nmod", t1);

            TIMEIT_START;
            _loop_redc_half_get_nmod(z, x, N, ctx_half);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_get_nmod", t1);
        }

        flint_printf("\n");

        TIMEIT_START;
        _loop_nmod_mul(z, x, y, N, mod);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_mul", t1);

        if (pbits == FLINT_BITS)
        {
            TIMEIT_START;
            _loop_nmod_mul_fullword(z, x, y, N, mod);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("_nmod_mul_fullword", t1);
        }

        TIMEIT_START;
        _loop_nmod_redc_mul(z, x_red, y_red, N, ctx);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redc_mul", t1);

        if (nmod_redc_can_use_fast(ctx))
        {
            TIMEIT_START;
            _loop_nmod_redc_fast_mul(z, x_red, y_red, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_fast_mul", t1);
        }

        if (pbits <= FLINT_BITS / 2 - 2)
        {
            TIMEIT_START;
            _loop_nmod_redc_half_mul(z, x_red_half, y_red_half, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_mul", t1);

            TIMEIT_START;
            _loop_nmod_redc_half_fast_mul(z, x_red_half_fast, y_red_half_fast, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_fast_mul", t1);
        }

        flint_printf("\n");

        TIMEIT_START;
        _loop_nmod_add(z, x, y, N, mod);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_add", t1);

        TIMEIT_START;
        _loop_nmod_redc_add(z, x_red, y_red, N, ctx);
        TIMEIT_STOP_VALUES(tcpu, t1);
        PRINTF2("nmod_redc_add", t1);

        if (nmod_redc_can_use_fast(ctx))
        {
            TIMEIT_START;
            _loop_nmod_redc_fast_add(z, x_red, y_red, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_fast_add", t1);

            TIMEIT_START;
            _loop_nmod_redc_fast_mul_two(z, x_red, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_fast_mul_two", t1);
        }

        if (pbits <= FLINT_BITS / 2 - 2)
        {
            TIMEIT_START;
            _loop_nmod_redc_half_add(z, x_red_half, y_red_half, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_add", t1);

            TIMEIT_START;
            _loop_nmod_redc_half_fast_add(z, x_red_half_fast, y_red_half_fast, N, ctx);
            TIMEIT_STOP_VALUES(tcpu, t1);
            PRINTF2("nmod_redc_half_fast_add", t1);
        }

        flint_printf("\n%10s  %24s  %24s  %24s  %24s  %24s  %24s\n", "ebits",
            "_nmod_pow_ui_binexp", "_nmod_redc_pow_ui", "_nmod_redc_fast_pow_ui",
            "_nmod_redc_half_pow_ui", "_nmod_redc_half_fast_pow_ui", "nmod_pow_ui");

        for (iebits = 0; (ebits = ebits_arr[iebits]) != 0; iebits++)
        {
            for (i = 0; i < N; i++)
            {
                exps[i] = n_randbits(state, ebits);
                exps[i] |= UWORD(1) << (ebits - 1);
            }

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += _nmod_pow_ui_binexp(x[i], exps[i], mod);
            }
            TIMEIT_STOP_VALUES(tcpu, t1);

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += _nmod_redc_pow_ui(x_red[i], exps[i], ctx);
            }
            TIMEIT_STOP_VALUES(tcpu, t2);

            if (nmod_redc_can_use_fast(ctx))
            {
                TIMEIT_START;
                for (i = 0; i < N; i++)
                {
                    s += _nmod_redc_fast_pow_ui(x_red[i], exps[i], ctx);
                }
                TIMEIT_STOP_VALUES(tcpu, t3);
            }
            else
                t3 = D_NAN;

            if (pbits <= FLINT_BITS / 2 - 2)
            {
                TIMEIT_START;
                for (i = 0; i < N; i++)
                {
                    s += _nmod_redc_half_pow_ui(x_red_half[i], exps[i], ctx);
                }
                TIMEIT_STOP_VALUES(tcpu, t4);

                TIMEIT_START;
                for (i = 0; i < N; i++)
                {
                    s += _nmod_redc_half_fast_pow_ui(x_red_half_fast[i], exps[i], ctx);
                }
                TIMEIT_STOP_VALUES(tcpu, t5);
            }
            else
            {
                t4 = t5 = D_NAN;
            }

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += nmod_pow_ui(x[i], exps[i], mod);
            }
            TIMEIT_STOP_VALUES(tcpu, t6);

            t1 /= N;
            t2 /= N;
            t3 /= N;
            t4 /= N;
            t5 /= N;
            t6 /= N;

            t1 /= 1e-9;
            t2 /= 1e-9;
            t3 /= 1e-9;
            t4 /= 1e-9;
            t5 /= 1e-9;
            t6 /= 1e-9;

            flint_printf("%10wu  %24g  %24g  %24g  %24g  %24g  %24g\n", ebits, t1, t2, t3, t4, t5, t6);
        }

        flint_printf("\n%10s  %24s  %24s  %24s  %24s\n", "ebits",
            "_nmod_2_pow_ui_binexp", "_nmod_redc_2_pow_ui", "_nmod_redc_fast_2_pow_ui", "nmod_2_pow_ui");

        for (iebits = 0; (ebits = ebits_arr[iebits]) != 0; iebits++)
        {
            if (ebits < 8)
                continue;

            for (i = 0; i < N; i++)
            {
                exps[i] = n_randbits(state, ebits);
                exps[i] |= UWORD(1) << (ebits - 1);
            }

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += _nmod_2_pow_ui_binexp(exps[i], mod);
            }
            TIMEIT_STOP_VALUES(tcpu, t1);

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += _nmod_redc_2_pow_ui(exps[i], ctx);
            }
            TIMEIT_STOP_VALUES(tcpu, t2);

            if (nmod_redc_can_use_fast(ctx))
            {
                TIMEIT_START;
                for (i = 0; i < N; i++)
                {
                    s += _nmod_redc_fast_2_pow_ui(exps[i], ctx);
                }
                TIMEIT_STOP_VALUES(tcpu, t3);
            }
            else
                t3 = D_NAN;

            TIMEIT_START;
            for (i = 0; i < N; i++)
            {
                s += nmod_2_pow_ui(exps[i], mod);
            }
            TIMEIT_STOP_VALUES(tcpu, t4);

            t1 /= N;
            t2 /= N;
            t3 /= N;
            t4 /= N;

            t1 /= 1e-9;
            t2 /= 1e-9;
            t3 /= 1e-9;
            t4 /= 1e-9;

            flint_printf("%10wu  %24g  %24g  %24g  %24g\n", ebits, t1, t2, t3, t4);
        }
    }

    /* Use value to prevent compiler from eliminating code */
    n_sqrt(s);

    _nmod_vec_clear(x);
    _nmod_vec_clear(y);
    _nmod_vec_clear(x_red);
    _nmod_vec_clear(y_red);
    _nmod_vec_clear(x_red_half);
    _nmod_vec_clear(y_red_half);
    _nmod_vec_clear(x_red_fast);
    _nmod_vec_clear(y_red_fast);
    _nmod_vec_clear(x_red_half_fast);
    _nmod_vec_clear(y_red_half_fast);
    _nmod_vec_clear(exps);

    flint_rand_clear(state);
}


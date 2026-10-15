/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "mpn_extras.h"
#include "nfloat.h"
#include "mp_real.h"
#include "gr.h"
#include "mp_real/impl.h"

/* The few-limb regimes touch no shared cache: with the tuned defaults,
   the elementary functions at up to MP_REAL_STATIC_REGIME_N limbs (and
   the nfloat functions at those precisions) and the constants pi/4,
   log 2, 2/pi, 1/log 2 - 1 up to MP_REAL_CONST_STATIC_N limbs read
   static tables only, so that a program working at these sizes builds
   no table at run time.  Up to MP_REAL_STATIC_TABLES_N limbs beyond
   that, the only run-time storage touched is the lazily filled dense
   tables (dense_rt.c, 64-bit only), which the guard limbs of the
   functions reach from n = MP_REAL_STATIC_REGIME_N + 1 (two guard
   limbs: expm1, sinh, cosh, tanh, asinh and the nfloat exp10, log10,
   log1p) and n = MP_REAL_STATIC_REGIME_N + 2 (one guard limb: the
   rest); the versioned caches -- the bitwise and Machin tables and the
   constants -- are never touched.  Arguments stay below 2^512 in
   magnitude: the trigonometric reduction of a larger argument
   legitimately needs pi beyond the static table. */

#if FLINT_BITS == 64
#define MP_REAL_STATIC_REGIME_N 8
#define MP_REAL_STATIC_TABLES_N 10
#else
/* the limb bitwise paths read the static prefixes of the bitwise tables
   (impl.h, r <= 32) at the tuned r = 32 of working precisions below 11
   (trigonometric) or 12 (exponential) limbs: at 9 limbs, the guard limbs
   of tanh, tan, sin_cos_pi and the nfloat trigonometric functions reach
   11 */
#define MP_REAL_STATIC_REGIME_N 8
#define MP_REAL_STATIC_TABLES_N 8
#endif

#if FLINT_BITS == 64
#define DENSE_RT_ALLOCATED() _mp_real_dense_rt_allocated()
#else
#define DENSE_RT_ALLOCATED() 0
#endif

static void
report(const char * what, slong n)
{
    int id;

    flint_printf("%s at n = %wd published a shared cache: exp_logs r = %wd, "
        "atans r = %wd, log_primes limbs = %wd, constants limbs:", what, n,
        _mp_real_exp_logs_max_index(), _mp_real_atans_max_index(),
        _mp_real_log_primes_max_limbs());
    for (id = 0; id < 12; id++)
        flint_printf(" %wd", _mp_real_const_cached_limbs(id));
#if FLINT_BITS == 64
    flint_printf(", dense_rt storages = %wd", DENSE_RT_ALLOCATED());
#endif
    flint_printf("\n");
}

#define CHECK(what) \
    do { \
        if (_flint_shared_cache_count() - DENSE_RT_ALLOCATED() != 0 \
            || (n <= MP_REAL_STATIC_REGIME_N && DENSE_RT_ALLOCATED() != 0)) \
        { \
            report(what, n); \
            fails++; \
        } \
        flint_cleanup_master(); \
    } while (0)

/* a random x with |x| < 2^emax (and |x| >= 2^emin unless zero) */
static void
random_mp_real(mp_real_t x, flint_rand_t state, slong n, slong emin, slong emax)
{
    ulong m[20];
    slong k, e;

    for (k = 0; k < n; k++)
        m[k] = n_randtest(state);
    m[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
    e = emin + n_randint(state, emax - emin + 1);
    _mp_real_set_mpn_2exp(x, m, n, e - FLINT_BITS * n);
    if (n_randint(state, 2))
        mp_real_neg(x, x);
}

TEST_FUNCTION_START(mp_real_static_regime, state)
{
    slong iter, n, k, fails = 0;

    /* every cache empty and unregistered */
    flint_cleanup_master();

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, y, z, w;
        slong prec;

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(z);
        mp_real_init(w);

        n = 1 + n_randint(state, MP_REAL_STATIC_TABLES_N);
        prec = FLINT_BITS * n - n_randint(state, FLINT_BITS);

        /* exponential family: |x| < 2^20 */
        random_mp_real(x, state, n, -FLINT_BITS * n - 20, 20);
        mp_real_exp_bits(y, x, prec);
        CHECK("exp");
        mp_real_expm1_bits(y, x, prec);
        CHECK("expm1");
        mp_real_sinh_cosh_bits(y, z, x, prec);
        CHECK("sinh_cosh");
        mp_real_tanh_bits(y, x, prec);
        CHECK("tanh");
        mp_real_asinh_bits(y, x, prec);
        CHECK("asinh");

        /* logarithms: x > 0, below 2^512 */
        random_mp_real(x, state, n, -512, 512);
        x->negative = 0;
        mp_real_log_bits(y, x, prec);
        CHECK("log");
        mp_real_log1p_bits(y, x, prec);
        CHECK("log1p");
        mp_real_sqrt(y, x, n);
        CHECK("sqrt");
        mp_real_root_ui(y, x, 3 + n_randint(state, 20), n);
        CHECK("root_ui");
        mp_real_acosh_bits(y, x, prec);
        CHECK("acosh");

        /* trigonometric: |x| < 2^512 */
        random_mp_real(x, state, n, -FLINT_BITS * n - 20, 512);
        mp_real_sin_cos_bits(y, z, x, prec);
        CHECK("sin_cos");
        mp_real_tan_bits(y, x, prec);
        CHECK("tan");
        mp_real_atan_bits(y, x, prec);
        CHECK("atan");
        mp_real_sin_cos_pi_bits(y, z, x, prec);
        CHECK("sin_cos_pi");
        mp_real_tan_pi_bits(y, x, prec);
        CHECK("tan_pi");
        mp_real_sin_cos_pi_ui_div_ui(y, z, n_randint(state, 100),
            1 + n_randint(state, 100), n);
        CHECK("sin_cos_pi_ui_div_ui");

        /* inverse trigonometric: |x| < 1 */
        random_mp_real(x, state, n, -FLINT_BITS * n - 20, 0);
        mp_real_asin_bits(y, x, prec);
        CHECK("asin");
        mp_real_acos_bits(y, x, prec);
        CHECK("acos");
        mp_real_atanh_bits(y, x, prec);
        CHECK("atanh");

        /* the constants with static tables, with the cache on */
        n = 1 + n_randint(state, MP_REAL_CONST_STATIC_N);
        mp_real_const_pi4(y, n, 1);
        CHECK("const_pi4");
        mp_real_const_log2(y, n, 1);
        CHECK("const_log2");
        mp_real_const_2_div_pi(y, n, 1);
        CHECK("const_2_div_pi");
        mp_real_const_inv_log2_frac(y, n, 1);
        CHECK("const_inv_log2_frac");

        /* nfloat at the same precisions */
        {
            gr_ctx_t ctx;
            nfloat_ptr a, b, c;
            double d;

            n = 1 + n_randint(state, MP_REAL_STATIC_TABLES_N);
            prec = FLINT_BITS * n;
            nfloat_ctx_init(ctx, prec, 0);
            a = gr_heap_init(ctx);
            b = gr_heap_init(ctx);
            c = gr_heap_init(ctx);

            d = (n_randint(state, 1000) + 1) / 7.0;
            if (n_randint(state, 2))
                d = -d;
            nfloat_set_d(a, d, ctx);
            k = (slong) n_randint(state, 80) - 70;
            nfloat_mul_2exp_si(a, a, k, ctx);

            nfloat_exp(b, a, ctx);
            CHECK("nfloat_exp");
            nfloat_expm1(b, a, ctx);
            CHECK("nfloat_expm1");
            nfloat_exp2(b, a, ctx);
            CHECK("nfloat_exp2");
            nfloat_exp10(b, a, ctx);
            CHECK("nfloat_exp10");
            nfloat_sin_cos(b, c, a, ctx);
            CHECK("nfloat_sin_cos");
            nfloat_tan(b, a, ctx);
            CHECK("nfloat_tan");
            nfloat_atan(b, a, ctx);
            CHECK("nfloat_atan");
            nfloat_sinh_cosh(b, c, a, ctx);
            CHECK("nfloat_sinh_cosh");
            nfloat_tanh(b, a, ctx);
            CHECK("nfloat_tanh");
            nfloat_sin_cos_pi(b, c, a, ctx);
            CHECK("nfloat_sin_cos_pi");
            nfloat_tan_pi(b, a, ctx);
            CHECK("nfloat_tan_pi");
            nfloat_abs(a, a, ctx);
            nfloat_log(b, a, ctx);
            CHECK("nfloat_log");
            nfloat_log1p(b, a, ctx);
            CHECK("nfloat_log1p");
            nfloat_log2(b, a, ctx);
            CHECK("nfloat_log2");
            nfloat_log10(b, a, ctx);
            CHECK("nfloat_log10");
            nfloat_sqrt(b, a, ctx);
            CHECK("nfloat_sqrt");
            nfloat_set_d(c, 1.25, ctx);
            nfloat_pow(b, a, c, ctx);
            CHECK("nfloat_pow");
            nfloat_atan2(b, a, c, ctx);
            CHECK("nfloat_atan2");
            nfloat_pi(b, ctx);
            CHECK("nfloat_pi");

            gr_heap_clear(a, ctx);
            gr_heap_clear(b, ctx);
            gr_heap_clear(c, ctx);
            gr_ctx_clear(ctx);
        }

        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(z);
        mp_real_clear(w);
    }

    if (fails)
        TEST_FUNCTION_FAIL("%wd evaluations touched a shared cache\n", fails);

    TEST_FUNCTION_END(state);
}

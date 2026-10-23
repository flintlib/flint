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
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

#if FLINT_BITS == 64

/* |(y, len) B^-n - c| <= err B^-n, c an arb_t (wide enough) */
static int
_within(nn_srcptr y, slong len, slong n, const arb_t c, ulong err)
{
    fmpz_t f;
    arb_t t;
    mag_t m;
    int ok;

    fmpz_init(f);
    arb_init(t);
    mag_init(m);
    fmpz_set_ui_array(f, y, len);
    arb_set_fmpz(t, f);
    arb_mul_2exp_si(t, t, -FLINT_BITS * n);
    arb_sub(t, t, c, FLINT_BITS * n + 128);
    arb_get_mag(m, t);
    mag_mul_2exp_si(m, m, FLINT_BITS * n);
    ok = (mag_cmp_2exp_si(m, 62) < 0) && (mag_get_d(m) <= (double) err);
    fmpz_clear(f);
    arb_clear(t);
    mag_clear(m);
    return ok;
}

/* inputs that hit the index boundaries of the seven-bit levels: a
   pattern of top bits (all equal within each level but one) over
   all-zero or all-one low limbs, or plain random limbs */
static void
_dense_input(nn_ptr x, slong n, flint_rand_t state, slong iter)
{
    slong i;

    switch (iter % 6)
    {
        case 0:
            flint_mpn_rrandom(x, state, n);
            break;
        case 1:
            for (i = 0; i < n; i++)
                x[i] = n_randlimb(state);
            break;
        case 2: case 3:
        {
            /* ten levels of seven bits, the last one straddling the
               top two limbs */
            ulong low = (iter % 6 == 2) ? 0 : ~UWORD(0), top = 0, g;
            slong k, b = n_randint(state, 10);

            for (k = 0; k < 10; k++)
            {
                g = (k == b) ? n_randint(state, 128)
                    : (n_randint(state, 2) ? 127 : 0);
                top = (k < 9) ? (top << 7) | g : (top << 1) | (g >> 6);
                if (k == 9)
                    low = (g << 58) | (low >> 6);
            }
            for (i = 0; i < n - 2; i++)
                x[i] = (iter % 6 == 2) ? 0 : ~UWORD(0);
            x[n - 2] = low;
            x[n - 1] = top;
            break;
        }
        case 4:
            /* sparse */
            for (i = 0; i < n; i++)
                x[i] = (n_randint(state, 3) == 0) ? n_randlimb(state)
                    >> n_randint(state, FLINT_BITS) : 0;
            break;
        default:
            for (i = 0; i < n; i++)
                x[i] = ~UWORD(0);
            if (n_randint(state, 2))
                x[n - 1] = n_randlimb(state) >> n_randint(state, 8);
            break;
    }
}

/* the run-time tables from several threads at once: each task evaluates
   all kernels at a size in the run-time range (the first calls race to
   build and widen the tables), against the same calls made serially */
#define RT_TASKS 48
#define RT_LEN (4 * MP_REAL_DENSE_RT_MAX + 4)

typedef struct
{
    ulong (* out)[RT_LEN];
}
_rt_args_t;

static slong
_rt_task_n(slong i)
{
    return MP_REAL_DENSE_W + 1 + (i * 7) % (MP_REAL_DENSE_RT_MAX - MP_REAL_DENSE_W);
}

static void
_rt_task(slong i, void * args)
{
    ulong * y = ((_rt_args_t *) args)->out[i];
    ulong x[MP_REAL_DENSE_RT_MAX], err;
    slong n = _rt_task_n(i), k;

    for (k = 0; k < n; k++)
        x[k] = UWORD(0x9e3779b97f4a7c15) * (ulong) (i * 31 + k + 1);
    flint_mpn_zero(y, RT_LEN);
    _mp_real_exp_opt(y, &err, x, n);
    _mp_real_log1p_opt(y + n + 1, &err, x, n);
    _mp_real_atan_opt(y + 2 * n + 1, &err, x, n);
    _mp_real_sin_cos_opt(y + 3 * n + 1, y + 2, &err, x, n);
}

#endif

TEST_FUNCTION_START(mp_real_dense_kernels, state)
{
#if FLINT_BITS == 64
    slong iter;

    if (MP_REAL_DENSE_RT_MAX > MP_REAL_DENSE_W)
    {
        ulong (* par)[RT_LEN] = flint_malloc(RT_TASKS * sizeof(*par));
        ulong (* ser)[RT_LEN] = flint_malloc(RT_TASKS * sizeof(*ser));
        _rt_args_t a;
        slong i, nthreads = flint_get_num_threads();

        flint_set_num_threads(4);
        a.out = par;
        flint_parallel_do(_rt_task, &a, RT_TASKS, 4, FLINT_PARALLEL_STRIDED);
        flint_set_num_threads(nthreads);
        a.out = ser;
        for (i = 0; i < RT_TASKS; i++)
            _rt_task(i, &a);
        for (i = 0; i < RT_TASKS; i++)
            if (mpn_cmp(par[i], ser[i], RT_LEN) != 0)
                TEST_FUNCTION_FAIL("threads: task %wd, n = %wd\n", i, _rt_task_n(i));

        flint_free(par);
        flint_free(ser);
    }

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong n = 4 + n_randint(state, MP_REAL_DENSE_RT_MAX - 3), prec = FLINT_BITS * n + 128;
        ulong x[MP_REAL_DENSE_RT_MAX + 1], y[MP_REAL_DENSE_RT_MAX + 3], y2[MP_REAL_DENSE_RT_MAX + 3], err;
        fmpz_t f;
        arb_t v, c, s;

        fmpz_init(f);
        arb_init(v);
        arb_init(c);
        arb_init(s);

        _dense_input(x, n, state, iter);
        fmpz_set_ui_array(f, x, n);
        arb_set_fmpz(v, f);
        arb_mul_2exp_si(v, v, -FLINT_BITS * n);

        if (!_mp_real_exp_opt(y, &err, x, n))
            TEST_FUNCTION_FAIL("exp: no kernel at n = %wd\n", n);
        arb_exp(c, v, prec);
        if (!_within(y, n + 1, n, c, err))
            TEST_FUNCTION_FAIL("exp: n = %wd, iter = %wd\n", n, iter);

        if (!_mp_real_log1p_opt(y, &err, x, n))
            TEST_FUNCTION_FAIL("log1p: no kernel at n = %wd\n", n);
        arb_log1p(c, v, prec);
        if (!_within(y, n, n, c, err))
            TEST_FUNCTION_FAIL("log1p: n = %wd, iter = %wd\n", n, iter);

        if (!_mp_real_atan_opt(y, &err, x, n))
            TEST_FUNCTION_FAIL("atan: no kernel at n = %wd\n", n);
        arb_atan(c, v, prec);
        if (!_within(y, n, n, c, err))
            TEST_FUNCTION_FAIL("atan: n = %wd, iter = %wd\n", n, iter);

        if (!_mp_real_sin_cos_opt(y, y2, &err, x, n))
            TEST_FUNCTION_FAIL("sin_cos: no kernel at n = %wd\n", n);
        arb_sin_cos(s, c, v, prec);
        if (!_within(y, n + 1, n, s, err) || !_within(y2, n + 1, n, c, err))
            TEST_FUNCTION_FAIL("sin_cos: n = %wd, iter = %wd\n", n, iter);

        /* the reciprocal: I <= 1/(1 + F), within six ulps */
        if (n <= MP_REAL_DENSE_W)
        {
            ulong I[MP_REAL_DENSE_W + 1];

            _flint_mpn_recip_small(I, x, n, FLINT_BITS * n - 3);
            arb_one(c);
            arb_add(c, c, v, prec);
            arb_inv(c, c, prec);
            fmpz_set_ui_array(f, I, n);
            arb_set_fmpz(s, f);
            arb_mul_2exp_si(s, s, -FLINT_BITS * n);
            arb_sub(s, c, s, prec);
            if (!arb_is_nonnegative(s) || !_within(I, n, n, c, 6))
            {
                slong i;
                flint_printf("F =");
                for (i = n - 1; i >= 0; i--)
                    flint_printf(" %016wx", x[i]);
                flint_printf("\nI =");
                for (i = n - 1; i >= 0; i--)
                    flint_printf(" %016wx", I[i]);
                flint_printf("\nd = "); arb_printd(s, 20); flint_printf("\n");
                TEST_FUNCTION_FAIL("recip: n = %wd, iter = %wd\n", n, iter);
            }
        }

        fmpz_clear(f);
        arb_clear(v);
        arb_clear(c);
        arb_clear(s);
    }
#endif

    TEST_FUNCTION_END(state);
}

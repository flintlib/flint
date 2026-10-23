/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mag.h"
#include "arf.h"
#include "decimal.h"
#include "fmpq.h"
#include "gr.h"

/* whether n = 10^k */
static int
fmpz_is_pow10(const fmpz_t n)
{
    fmpz_t t;
    int res;
    fmpz_init_set(t, n);
    while (fmpz_divisible_ui(t, 10) && !fmpz_is_zero(t))
        fmpz_divexact_ui(t, t, 10);
    res = fmpz_is_one(t);
    fmpz_clear(t);
    return res;
}

/* random finite decmag with a small exponent */
static void
_decmag_randtest_small(decmag_t x, flint_rand_t state, gr_ctx_t ctx)
{
    ulong m = 1 + n_randint(state, DECIMAL_CTX_RAD_POW(ctx) - 1);
    slong e = (slong) n_randint(state, 61) - 30;

    if (n_randint(state, 8) == 0)
        _decmag_zero(x, ctx);
    else
        _decmag_set_ui_10exp_si(x, m, e, ctx);
}

/* checks lo <= exact <= hi and that hi <= exact * (1 + slack * 10^(1-rp))
   (tightness) when exact is nonzero */
static void
_check_bounds(const char * name, const decmag_t lo, const decmag_t hi, const fmpq_t exact,
    slong slack, gr_ctx_t ctx)
{
    fmpq_t l, h, t;

    fmpq_init(l);
    fmpq_init(h);
    fmpq_init(t);

    _decmag_get_fmpq(l, lo, ctx);
    _decmag_get_fmpq(h, hi, ctx);

    if (fmpq_cmp(l, exact) > 0 || fmpq_cmp(h, exact) < 0)
    {
        flint_printf("FAIL: %s bounds\n", name);
        flint_printf("lo = %{fmpq}\nexact = %{fmpq}\nhi = %{fmpq}\n", l, exact, h);
        flint_abort();
    }

    /* tightness: hi <= exact * (1 + slack * 10^(1-rp)) and lo >= exact * (1 - slack * 10^(1-rp)) */
    if (!fmpq_is_zero(exact))
    {
        fmpq_t eps;
        fmpq_init(eps);
        fmpq_set_si(eps, slack, 1);
        fmpz_ui_pow_ui(fmpq_denref(eps), 10, DECIMAL_CTX_RAD_PREC(ctx) - 1);
        fmpq_canonicalise(eps);

        fmpq_add_si(t, eps, 1);
        fmpq_mul(t, t, exact);
        if (fmpq_cmp(h, t) > 0)
        {
            flint_printf("FAIL: %s upper bound not tight\n", name);
            flint_printf("exact = %{fmpq}\nhi = %{fmpq}\n", exact, h);
            flint_abort();
        }

        fmpq_sub_si(t, eps, 1);
        fmpq_neg(t, t);
        fmpq_mul(t, t, exact);
        if (!_decmag_is_zero(lo, ctx) && fmpq_cmp(l, t) < 0)
        {
            flint_printf("FAIL: %s lower bound not tight\n", name);
            flint_printf("exact = %{fmpq}\nlo = %{fmpq}\n", exact, l);
            flint_abort();
        }

        fmpq_clear(eps);
    }

    fmpq_clear(l);
    fmpq_clear(h);
    fmpq_clear(t);
}

TEST_FUNCTION_START(decmag, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, ctx0;
        decmag_t x, y, lo, hi;
        fmpq_t qx, qy, qz;
        decfloat_t f;
        slong slack;

        gr_ctx_init_decball_randtest(ctx, state, 30);
        /* the operands may come from a context with another radius precision,
           in which case the results are slightly less tight */
        gr_ctx_init_decball_randtest(ctx0, state, 30);
        if (n_randint(state, 2))
            decimal_ctx_set_rad_prec(ctx0, DECIMAL_CTX_RAD_PREC(ctx));
        slack = (DECIMAL_CTX_RAD_PREC(ctx0) == DECIMAL_CTX_RAD_PREC(ctx)) ? 2 : 4;

        _decmag_init(x, ctx);
        _decmag_init(y, ctx);
        _decmag_init(lo, ctx);
        _decmag_init(hi, ctx);
        fmpq_init(qx);
        fmpq_init(qy);
        fmpq_init(qz);
        decfloat_init(f, ctx);

        _decmag_randtest_small(x, state, ctx0);
        _decmag_randtest_small(y, state, ctx0);

        _decmag_get_fmpq(qx, x, ctx);
        _decmag_get_fmpq(qy, y, ctx);

        /* add */
        fmpq_add(qz, qx, qy);
        _decmag_add(hi, x, y, ctx);
        _decmag_add_lower(lo, x, y, ctx);
        _check_bounds("add", lo, hi, qz, slack, ctx);

        /* aliasing */
        _decmag_set(lo, x, ctx);
        _decmag_add(lo, lo, y, ctx);
        if (!_decmag_equal(lo, hi, ctx)) { flint_printf("FAIL: add aliasing\n"); flint_abort(); }
        _decmag_set(lo, y, ctx);
        _decmag_add(lo, x, lo, ctx);
        if (!_decmag_equal(lo, hi, ctx)) { flint_printf("FAIL: add aliasing 2\n"); flint_abort(); }

        /* sub_lower */
        fmpq_sub(qz, qx, qy);
        if (fmpq_sgn(qz) < 0)
            fmpq_zero(qz);
        _decmag_sub_lower(lo, x, y, ctx);
        {
            fmpq_t l;
            fmpq_init(l);
            _decmag_get_fmpq(l, lo, ctx);
            if (fmpq_cmp(l, qz) > 0)
            {
                flint_printf("FAIL: sub_lower\n");
                flint_printf("x = %{fmpq}, y = %{fmpq}, lo = %{fmpq}\n", qx, qy, l);
                flint_abort();
            }
            /* tightness when x > y significantly */
            if (!fmpq_is_zero(qz))
            {
                fmpq_t t, eps;
                fmpq_init(t);
                fmpq_init(eps);
                fmpq_set_si(eps, slack, 1);
                fmpz_ui_pow_ui(fmpq_denref(eps), 10, DECIMAL_CTX_RAD_PREC(ctx) - 1);
                fmpq_canonicalise(eps);
                /* lo >= (x - y) - eps * x */
                fmpq_mul(t, eps, qx);
                fmpq_sub(t, qz, t);
                if (fmpq_cmp(l, t) < 0)
                {
                    flint_printf("FAIL: sub_lower not tight\n");
                    flint_printf("x = %{fmpq}, y = %{fmpq}, lo = %{fmpq}\n", qx, qy, l);
                    flint_abort();
                }
                fmpq_clear(t);
                fmpq_clear(eps);
            }
            fmpq_clear(l);
        }

        /* mul */
        fmpq_mul(qz, qx, qy);
        _decmag_mul(hi, x, y, ctx);
        _decmag_mul_lower(lo, x, y, ctx);
        _check_bounds("mul", lo, hi, qz, slack, ctx);

        /* div */
        if (!fmpq_is_zero(qy))
        {
            fmpq_div(qz, qx, qy);
            _decmag_div(hi, x, y, ctx);
            _decmag_div_lower(lo, x, y, ctx);
            _check_bounds("div", lo, hi, qz, slack, ctx);
        }

        /* sqrt: check lo^2 <= x <= hi^2 */
        {
            fmpq_t l, h;
            fmpq_init(l);
            fmpq_init(h);
            _decmag_sqrt(hi, x, ctx);
            _decmag_sqrt_lower(lo, x, ctx);
            _decmag_get_fmpq(l, lo, ctx);
            _decmag_get_fmpq(h, hi, ctx);
            fmpq_mul(l, l, l);
            fmpq_mul(h, h, h);
            if (fmpq_cmp(l, qx) > 0 || fmpq_cmp(h, qx) < 0)
            {
                flint_printf("FAIL: sqrt bounds\n");
                flint_printf("x = %{fmpq}\nlo^2 = %{fmpq}\nhi^2 = %{fmpq}\n", qx, l, h);
                flint_abort();
            }
            fmpq_clear(l);
            fmpq_clear(h);
        }

        /* conversion from decfloat: lower <= |f| <= upper */
        GR_MUST_SUCCEED(decfloat_randtest(f, state, ctx));
        if (decfloat_get_fmpq(qz, f, ctx) == GR_SUCCESS)
        {
            fmpq_abs(qz, qz);
            _decmag_set_decfloat(hi, f, ctx);
            _decmag_set_decfloat_lower(lo, f, ctx);
            _check_bounds("set_decfloat", lo, hi, qz, 2, ctx);

            /* ulp */
            if (!fmpq_is_zero(qz))
            {
                slong prec = 1 + n_randint(state, 30);
                slong E;
                fmpq_t u, v;
                fmpq_init(u);
                fmpq_init(v);
                _decmag_set_ulp(hi, f, prec, ctx);
                _decmag_get_fmpq(u, hi, ctx);
                GR_MUST_SUCCEED(decfloat_get_sci_exp_si(&E, f, ctx) ? GR_SUCCESS : GR_UNABLE);
                /* u == 10^(E - prec + 1) */
                fmpq_one(v);
                if (E - prec + 1 >= 0)
                    fmpz_ui_pow_ui(fmpq_numref(v), 10, E - prec + 1);
                else
                    fmpz_ui_pow_ui(fmpq_denref(v), 10, prec - 1 - E);
                if (!fmpq_equal(u, v))
                {
                    flint_printf("FAIL: ulp\n");
                    flint_printf("f = %{gr}, prec = %wd, u = %{fmpq}\n", f, ctx, prec, u);
                    flint_abort();
                }
                fmpq_clear(u);
                fmpq_clear(v);
            }
        }

        /* exact conversion to decfloat and back */
        _decmag_set(lo, x, ctx);
        GR_MUST_SUCCEED(_decmag_get_decfloat(f, x, ctx));
        _decmag_set_decfloat(hi, f, ctx);
        if (!_decmag_equal(lo, hi, ctx))
        {
            flint_printf("FAIL: decfloat roundtrip\n");
            flint_printf("x = %{fmpq}, f = %{gr}\n", qx, f, ctx);
            flint_abort();
        }

        /* cmp, max, min, is_10exp: values, independent of the mantissa digits */
        {
            int c1 = _decmag_cmp(x, y, ctx);
            int c2 = fmpq_cmp(qx, qy);
            /* max rounds up and min rounds down to the radius precision */
            _decmag_max(hi, x, y, ctx);
            _decmag_min(lo, x, y, ctx);
            if (c1 != c2 || _decmag_cmp(hi, x, ctx) < 0 || _decmag_cmp(hi, y, ctx) < 0
                || _decmag_cmp(lo, x, ctx) > 0 || _decmag_cmp(lo, y, ctx) > 0
                || (slack == 2 && !_decmag_equal(hi, (c1 >= 0) ? x : y, ctx))
                || (slack == 2 && !_decmag_equal(lo, (c1 <= 0) ? x : y, ctx))
                || (_decmag_equal(x, y, ctx) != (c1 == 0)))
            {
                flint_printf("FAIL: cmp\n");
                flint_printf("x = %{fmpq} (m = %wu, e = %{fmpz}), y = %{fmpq} (m = %wu, e = %{fmpz}), c1 = %d, c2 = %d\n", qx, x->m, &x->exp, qy, y->m, &y->exp, c1, c2);
                gr_ctx_println(ctx);
                flint_abort();
            }
        }

        if (!_decmag_is_zero(x, ctx))
        {
            int c1 = _decmag_is_10exp(x, ctx);
            int c2 = fmpz_is_one(fmpq_numref(qx)) ? fmpz_is_pow10(fmpq_denref(qx)) : (fmpz_is_one(fmpq_denref(qx)) && fmpz_is_pow10(fmpq_numref(qx)));
            if (c1 != c2)
            {
                flint_printf("FAIL: is_10exp\n");
                flint_printf("x = %{fmpq}, c1 = %d, c2 = %d\n", qx, c1, c2);
                flint_abort();
            }
        }

        _decmag_clear(x, ctx);
        _decmag_clear(y, ctx);
        _decmag_clear(lo, ctx);
        _decmag_clear(hi, ctx);
        fmpq_clear(qx);
        fmpq_clear(qy);
        fmpq_clear(qz);
        decfloat_clear(f, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(ctx0);
    }

    /* set_mag: upper bound, tight to the radius precision */
    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decmag_t x;
        mag_t m;
        arf_t a, b;
        decfloat_t f;
        slong prec_bits = 200;
        int which = n_randint(state, 4);

        gr_ctx_init_decball_randtest(ctx, state, 30);
        _decmag_init(x, ctx);
        mag_init(m);
        arf_init(a);
        arf_init(b);
        decfloat_init(f, ctx);

        mag_randtest_special(m, state, which == 0 ? 10 : (which == 1 ? 100 : (which == 2 ? 1000 : 3000)));

        _decmag_set_mag(x, m, ctx);

        if (mag_is_zero(m) || mag_is_inf(m))
        {
            if (mag_is_zero(m) != DECMAG_IS_ZERO(x) || mag_is_inf(m) != DECMAG_IS_INF(x))
            {
                flint_printf("FAIL: set_mag special\n");
                flint_abort();
            }
        }
        else
        {
            /* mag(x) >= m, mag(x) <= m (1 + 10^(1-rad_prec)) (1 + 2^-30) */
            arf_set_mag(a, m);
            GR_MUST_SUCCEED(_decmag_get_decfloat(f, x, ctx));
            GR_MUST_SUCCEED(decfloat_get_arf(b, f, prec_bits, ARF_RND_DOWN, ctx));
            if (arf_cmp(b, a) < 0)
            {
                flint_printf("FAIL: set_mag lower\n");
                flint_printf("m = %{mag}\nx = %s\n", m, _decmag_get_str(x, ctx));
                flint_abort();
            }
            GR_MUST_SUCCEED(decfloat_get_arf(b, f, prec_bits, ARF_RND_UP, ctx));
            arf_mul_ui(a, a, 1000000000, prec_bits, ARF_RND_UP);
            arf_add_ui(a, a, 0, prec_bits, ARF_RND_UP);
            {
                /* a = m (10^9 + 10^(10 - rad_prec) + 1) 10^-9 */
                arf_t c;
                arf_init(c);
                arf_set_mag(c, m);
                arf_mul_ui(c, c, n_pow(10, 10 - DECIMAL_CTX_RAD_PREC(ctx)) + 1, prec_bits, ARF_RND_UP);
                arf_add(a, a, c, prec_bits, ARF_RND_UP);
                arf_div_ui(a, a, 1000000000, prec_bits, ARF_RND_UP);
                arf_clear(c);
            }
            if (arf_cmp(b, a) > 0)
            {
                flint_printf("FAIL: set_mag upper\n");
                flint_printf("m = %{mag}\nx = %s\nrad_prec = %wd\n", m, _decmag_get_str(x, ctx), DECIMAL_CTX_RAD_PREC(ctx));
                flint_abort();
            }
        }

        _decmag_clear(x, ctx);
        mag_clear(m);
        arf_clear(a);
        arf_clear(b);
        decfloat_clear(f, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

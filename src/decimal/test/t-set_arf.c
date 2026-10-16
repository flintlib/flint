/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "arb.h"
#include "decimal.h"
#include "gr.h"

/* random arf with a mantissa of up to mbits bits and an exponent of
   up to ebits bits, occasionally an exact decimal number */
static void
_randtest_arf(arf_t f, flint_rand_t state, slong mbits, slong ebits)
{
    if (n_randint(state, 4) == 0)
    {
        /* d 10^k */
        fmpz_t d, t;
        fmpz_init(d);
        fmpz_init(t);
        fmpz_randtest_not_zero(d, state, 1 + n_randint(state, 60));
        fmpz_ui_pow_ui(t, 10, n_randint(state, 5000));
        fmpz_mul(d, d, t);
        arf_set_fmpz(f, d);
        if (n_randint(state, 2))
            arf_mul_2exp_si(f, f, -(slong) n_randint(state, 100));
        fmpz_clear(d);
        fmpz_clear(t);
    }
    else
    {
        arf_randtest_not_zero(f, state, 1 + n_randint(state, mbits), 1 + n_randint(state, ebits));
    }
}

TEST_FUNCTION_START(decfloat_set_arf, state)
{
    slong iter;

    /* correctly rounded conversion vs exact conversion followed by rounding */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t x, y, z, d;
        decmag_t err, err2;
        decimal_rounding_info info;
        arf_t f;
        fmpz_t m, t;
        slong prec;
        int rnd, status, status2;

        gr_ctx_init_decfloat_randtest(ctx, state, 60);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        decimal_ctx_set_rad_prec(ctx, 1 + n_randint(state, 9));
        prec = 1 + n_randint(state, 200);
        rnd = n_randint(state, 7);

        decfloat_init(x, ctx);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);
        decfloat_init(d, ctx);
        _decmag_init(err, ctx);
        _decmag_init(err2, ctx);
        arf_init(f);
        fmpz_init(m);
        fmpz_init(t);

        _randtest_arf(f, state, 2000, 15);
        arf_get_fmpz_2exp(m, t, f);

        status = _decfloat_set_round_fmpz_2exp_err(x, m, t, prec, rnd, &info, err, ctx);
        status2 = decfloat_set_round_arf(z, f, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);

        if (status != GR_SUCCESS || status2 != GR_SUCCESS)
        {
            flint_printf("FAIL: status %d %d\n", status, status2);
            flint_printf("f = %{arf}\n", f);
            flint_abort();
        }

        GR_MUST_SUCCEED(decfloat_set_round_info(y, z, prec, rnd, NULL, err2, ctx));

        if (decfloat_equal(x, y, ctx) != T_TRUE || decfloat_digits(x, ctx) > prec)
        {
            flint_printf("FAIL: rounding\n");
            flint_printf("f = %{arf}\nprec = %wd rnd = %d\nx = %{gr}\ny = %{gr}\nz = %{gr}\n", f, prec, rnd, x, ctx, y, ctx, z, ctx);
            flint_abort();
        }

        /* exact error e = |z - x| satisfies e <= err <= e (1 + 10^-(rad_prec-1)) + tiny */
        GR_MUST_SUCCEED(_decfloat_add(d, z, x, 1, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx));

        if (info.inexact != (decfloat_equal(x, z, ctx) != T_TRUE))
        {
            flint_printf("FAIL: inexact flag\n");
            flint_printf("f = %{arf}\nx = %{gr}\nz = %{gr}\n", f, x, ctx, z, ctx);
            flint_abort();
        }

        {
            decmag_t lo, hi;
            _decmag_init(lo, ctx);
            _decmag_init(hi, ctx);
            _decmag_set_decfloat_lower(lo, d, ctx);
            _decmag_set_decfloat(hi, d, ctx);
            /* hi = d (1 + 2 10^(1 - r)) plus the sticky slack below the guard
               digits, where r = min(rp, e + 1): the error is rounded up to rp
               digits, but only determined from the two leading limbs of the
               discarded part, the first of which may have a single digit
               (with upward rounding throughout) */
            {
                slong r = FLINT_MIN(DECIMAL_CTX_RAD_PREC(ctx), DECIMAL_CTX_E(ctx) + 1);
                _decmag_mul_ui(hi, hi, n_pow(10, r - 1) + 2, ctx);
                _decmag_div_ui(hi, hi, n_pow(10, r - 1), ctx);
            }
            {
                decmag_t u;
                _decmag_init(u, ctx);
                _decmag_set_ulp(u, x, prec + 19, ctx);    /* 10^(E - prec - 18) */
                _decmag_add(hi, hi, u, ctx);
                _decmag_clear(u, ctx);
            }
            if (_decmag_cmp(err, lo, ctx) < 0 || (DECIMAL_CTX_RAD_PREC(ctx) >= 2 && _decmag_cmp(err, hi, ctx) > 0))
            {
                flint_printf("FAIL: error bound\n");
                flint_printf("f = %{arf}\nprec = %wd rnd = %d rad_prec = %wd\nx = %{gr}\nz = %{gr}\nd = %{gr}\nerr = %s\nlo = %s\nhi = %s\n", f, prec, rnd, DECIMAL_CTX_RAD_PREC(ctx), x, ctx, z, ctx, d, ctx, _decmag_get_str(err, ctx), _decmag_get_str(lo, ctx), _decmag_get_str(hi, ctx));
                flint_abort();
            }
            _decmag_clear(lo, ctx);
            _decmag_clear(hi, ctx);
        }

        decfloat_clear(x, ctx);
        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        decfloat_clear(d, ctx);
        _decmag_clear(err, ctx);
        _decmag_clear(err2, ctx);
        arf_clear(f);
        fmpz_clear(m);
        fmpz_clear(t);
        gr_ctx_clear(ctx);
    }

    /* huge exponents: exact decimal detection and containment */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, bctx;
        decfloat_t x;
        decball_t X;
        arf_t f;
        arb_t a, b;
        fmpz_t d, t;
        slong prec, k;
        int exact, status;

        gr_ctx_init_decfloat_randtest(ctx, state, 30);
        decimal_ctx_set_exp_limits(ctx, WORD_MIN, WORD_MAX);
        prec = DECIMAL_CTX_PREC(ctx);
        if (prec == DECIMAL_PREC_EXACT)
        {
            prec = 1 + n_randint(state, 30);
            decimal_ctx_set_prec(ctx, prec);
        }
        _gr_ctx_init_decimal(bctx, DECIMAL_CTX_BALL, DECIMAL_CTX_E(ctx), prec, DECIMAL_CTX_RND(ctx), DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);

        decfloat_init(x, ctx);
        decball_init(X, bctx);
        arf_init(f);
        arb_init(a);
        arb_init(b);
        fmpz_init(d);
        fmpz_init(t);

        exact = n_randint(state, 2);

        if (exact)
        {
            /* d 10^k with d of at most prec digits (sometimes more) */
            slong nd = 1 + n_randint(state, prec + 2);
            fmpz_randtest_not_zero(d, state, (slong) (nd * 3.32) + 1);
            k = 1000 + n_randint(state, 30000);
            fmpz_ui_pow_ui(t, 10, k);
            fmpz_mul(d, d, t);
            arf_set_fmpz(f, d);
        }
        else
        {
            arf_randtest_not_zero(f, state, 1 + n_randint(state, 200), 1 + n_randint(state, 40));
        }

        status = decfloat_set_round_arf(x, f, prec, DECIMAL_CTX_RND(ctx), ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: huge conversion status %d\n", status);
            flint_printf("f = %{arf}\n", f);
            flint_abort();
        }

        /* the result is within one ulp of the true value */
        {
            arb_t u;
            mag_t md, mu;
            decmag_t ulp;
            arb_init(u);
            mag_init(md);
            mag_init(mu);
            _decmag_init(ulp, ctx);
            arb_set_arf(u, f);
            GR_MUST_SUCCEED(decfloat_get_arb(a, x, 3 * prec + 100, ctx));
            arb_sub(a, a, u, 3 * prec + 100);
            arb_get_mag(md, a);
            _decmag_set_ulp(ulp, x, prec, ctx);
            _decmag_get_mag(mu, ulp, ctx);
            if (mag_cmp(md, mu) > 0)
            {
                flint_printf("FAIL: huge conversion\n");
                flint_printf("f = %{arf}\nx = %{gr}\n", f, x, ctx);
                flint_abort();
            }
            arb_clear(u);
            mag_clear(md);
            mag_clear(mu);
            _decmag_clear(ulp, ctx);
        }

        if (exact)
        {
            decfloat_t z;
            decfloat_init(z, ctx);
            /* d has at most prec digits iff the exact conversion has */
            if (decfloat_set_round_arf(z, f, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx) == GR_SUCCESS)
            {
                if (decfloat_digits(z, ctx) <= prec && decfloat_equal(x, z, ctx) != T_TRUE)
                {
                    flint_printf("FAIL: exact huge conversion\n");
                    flint_printf("f = %{arf}\nx = %{gr}\nz = %{gr}\n", f, x, ctx, z, ctx);
                    flint_abort();
                }
            }
            decfloat_clear(z, ctx);
        }

        /* ball conversion: containment and tightness */
        arb_set_arf(a, f);
        if (n_randint(state, 2))
            arb_add_error_2exp_si(a, arf_abs_bound_lt_2exp_si(f) - 1 - n_randint(state, 500));

        status = decball_set_arb(X, a, bctx);
        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: set_arb status %d\n", status);
            flint_printf("a = %{arb}\n", a);
            flint_abort();
        }

        GR_MUST_SUCCEED(decball_get_arb(b, X, 3 * prec + 200, bctx));
        if (!arb_contains(b, a))
        {
            flint_printf("FAIL: set_arb containment\n");
            flint_printf("a = %{arb}\nX = %{gr}\nb = %{arb}\n", a, X, bctx, b);
            flint_abort();
        }

        if (arb_is_exact(a) && exact)
        {
            decfloat_t z;
            decfloat_init(z, ctx);
            if (decfloat_set_round_arf(z, f, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx) == GR_SUCCESS
                && decfloat_digits(z, ctx) <= prec && !_decball_is_exact(X, bctx))
            {
                flint_printf("FAIL: set_arb exactness\n");
                flint_printf("a = %{arb}\nX = %{gr}\n", a, X, bctx);
                flint_abort();
            }
            decfloat_clear(z, ctx);
        }

        /* the radius is at most the input radius plus one ulp, with a tiny slack */
        {
            decmag_t r, u;
            _decmag_init(r, bctx);
            _decmag_init(u, bctx);
            _decmag_set_mag(r, arb_radref(a), bctx);
            _decmag_set_ulp(u, &X->mid, prec, bctx);
            _decmag_add(r, r, u, bctx);
            _decmag_mul_ui(r, r, 1001, bctx);
            _decmag_div_ui(r, r, 1000, bctx);
            if (_decmag_cmp(&X->rad, r, bctx) > 0)
            {
                flint_printf("FAIL: set_arb radius\n");
                flint_printf("a = %{arb}\nX = %{gr}\nr = %s\n", a, X, bctx, _decmag_get_str(r, bctx));
                flint_abort();
            }
            _decmag_clear(r, bctx);
            _decmag_clear(u, bctx);
        }

        decfloat_clear(x, ctx);
        decball_clear(X, bctx);
        arf_clear(f);
        arb_clear(a);
        arb_clear(b);
        fmpz_clear(d);
        fmpz_clear(t);
        gr_ctx_clear(ctx);
        gr_ctx_clear(bctx);
    }

    TEST_FUNCTION_END(state);
}

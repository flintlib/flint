/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Spot checks of C-level helper functions that are not reachable through
   the generic gr interface (and hence not from the Python tests). */

#include <string.h>
#include <math.h>
#include "test_helpers.h"
#include "decimal.h"
#include "fmpq.h"
#include "arb.h"
#include "gr.h"
#include "double_extras.h"

static void
_check(int ok, const char * what)
{
    if (!ok)
    {
        flint_printf("FAIL: %s\n", what);
        flint_abort();
    }
}

static void
_check_str(char * s, const char * expected, const char * what)
{
    if (strcmp(s, expected))
    {
        flint_printf("FAIL: %s\ngot %s, expected %s\n", what, s, expected);
        flint_abort();
    }
    flint_free(s);
}

#define CHECK(expr) _check((expr), #expr)
#define CHECK_STR(s, t) _check_str((s), (t), #s)

static int
_arb_exp(arb_t res, const arb_t x, slong prec)
{
    arb_exp(res, x, prec);
    return GR_SUCCESS;
}

TEST_FUNCTION_START(decimal_api, state)
{
    gr_ctx_t F, B, C, CB;
    decmag_t u, v;
    decfloat_t x, y;
    decball_t X, Y;
    deccfloat_t z, w;
    deccball_t Z, W;
    fmpz_t n;
    fmpq_t q;
    gr_stream_t out;
    double d;
    slong i;

    gr_ctx_init_decfloat(F, 10, 0);
    gr_ctx_init_decball(B, 10, 0);
    gr_ctx_init_deccfloat(C, 10, 0);
    gr_ctx_init_deccball(CB, 10, 0);

    _decmag_init(u, B);
    _decmag_init(v, B);
    decfloat_init(x, F);
    decfloat_init(y, F);
    decball_init(X, B);
    decball_init(Y, B);
    deccfloat_init(z, C);
    deccfloat_init(w, C);
    deccball_init(Z, CB);
    deccball_init(W, CB);
    fmpz_init(n);
    fmpq_init(q);

    /* decmag: upper bounds, with rad_prec significant digits */
    _decmag_set_ui(u, 3, B);
    _decmag_inv(v, u, B);
    CHECK_STR(_decmag_get_str(v, B), "0.3334");
    _decmag_rsqrt(v, u, B);
    CHECK_STR(_decmag_get_str(v, B), "0.5774");
    _decmag_mul_10exp_si(v, u, -7, B);
    CHECK_STR(_decmag_get_str(v, B), "3e-7");
    CHECK(_decmag_cmp_10exp_si(v, -7, B) > 0 && _decmag_cmp_10exp_si(v, -6, B) < 0);
    CHECK(!_decmag_is_10exp(v, B));
    _decmag_min(v, u, v, B);
    CHECK_STR(_decmag_get_str(v, B), "3e-7");
    _decmag_mul_10exp_si(v, v, 7, B);
    CHECK(_decmag_equal(u, v, B));
    _decmag_set_10exp_si(v, 5, B);
    CHECK(_decmag_is_10exp(v, B) && _decmag_cmp_10exp_si(v, 5, B) == 0);

    fmpz_set_str(n, "123456789", 10);
    _decmag_set_fmpz(u, n, B);
    CHECK_STR(_decmag_get_str(u, B), "1.235e8");
    _decmag_set_fmpz_lower(u, n, B);
    CHECK_STR(_decmag_get_str(u, B), "1.234e8");

    CHECK(_decmag_set_d(u, 0.1, B) == GR_SUCCESS);
    CHECK_STR(_decmag_get_str(u, B), "0.1001");
    CHECK(_decmag_get_d(&d, u, B) == GR_SUCCESS && fabs(d - 0.1001) < 1e-15);
    CHECK(_decmag_set_d(v, -0.1, B) == GR_SUCCESS && _decmag_equal(u, v, B) && _decmag_set_d(v, D_NAN, B) == GR_DOMAIN);

    gr_stream_init_str(out);
    CHECK(_decmag_write(out, u, B) == GR_SUCCESS);
    CHECK_STR(out->s, "0.1001");

    for (i = 0; i < 100; i++)
    {
        _decmag_randtest(u, state, B);
        _decmag_randtest_special(v, state, B);
        _decmag_max(v, u, v, B);
        CHECK(_decmag_cmp(u, v, B) <= 0);
    }

    /* decmag: large integers and special values */
    fmpz_set_str(n, "1000000000000000000000000000001", 10);
    _decmag_set_fmpz(u, n, B);
    CHECK_STR(_decmag_get_str(u, B), "1.001e30");
    _decmag_set_fmpz_lower(u, n, B);
    CHECK_STR(_decmag_get_str(u, B), "1e30");
    fmpz_zero(n);
    _decmag_set_fmpz(u, n, B);
    _decmag_set_fmpz_lower(v, n, B);
    CHECK(_decmag_is_zero(u, B) && _decmag_is_zero(v, B));
    CHECK_STR(_decmag_get_str(u, B), "0");
    _decmag_mul_10exp_si(v, u, 5, B);
    _decmag_set_ui_10exp_fmpz_lower(v, 0, n, B);
    CHECK(_decmag_is_zero(v, B));
    _decmag_pow_ui(v, u, 0, B);
    _decmag_min(u, u, v, B);
    CHECK(_decmag_is_zero(u, B) && _decmag_cmp_10exp_si(v, 0, B) == 0);
    _decmag_div(v, v, u, B);
    CHECK(_decmag_is_inf(v, B));
    _decmag_div(u, u, v, B);
    _decmag_sub_lower(u, v, u, B);
    CHECK(_decmag_is_inf(u, B));
    CHECK(_decmag_set_d(v, D_INF, B) == GR_SUCCESS && _decmag_is_inf(v, B));
    /* an infinite radius converts to a decfloat only where infinities exist */
    CHECK(_decmag_get_decfloat(x, v, B) == GR_DOMAIN);
    {
        gr_ctx_t Fi;
        gr_ctx_init_decfloat(Fi, 10, DECIMAL_ALLOW_INF);
        CHECK(_decmag_get_decfloat(x, v, Fi) == GR_SUCCESS && DECFLOAT_IS_POS_INF(x));
        _decmag_set_decfloat_lower(u, x, B);
        CHECK(_decmag_is_inf(u, B));
        gr_ctx_clear(Fi);
    }

    /* decfloat */
    CHECK(decfloat_set_si_10exp_si(x, -25, -1, F) == GR_SUCCESS);
    CHECK_STR(decfloat_get_str(x, F), "-2.5");
    CHECK(_decfloat_cmp_ui(x, 2, F) < 0 && _decfloat_cmpabs_ui(x, 2, F) > 0 && _decfloat_cmpabs_ui(x, 3, F) < 0);
    CHECK(decfloat_get_fmpz_fixed_si(n, x, 0, DECIMAL_RND_FLOOR, F) == GR_SUCCESS && fmpz_equal_si(n, -3));
    CHECK(decfloat_get_fmpz_fixed_si(n, x, -3, DECIMAL_RND_NEAR, F) == GR_SUCCESS && fmpz_equal_si(n, -2500));
    CHECK(decfloat_get_fmpz_fixed_si(n, x, 1, DECIMAL_RND_UP, F) == GR_SUCCESS && fmpz_equal_si(n, -1));
    CHECK(decfloat_mul_10exp_si(y, x, 100, F) == GR_SUCCESS);
    CHECK_STR(decfloat_get_str(y, F), "-2.5e100");
    fmpz_set_si(n, 7);
    CHECK(decfloat_set_fmpz_10exp_si(y, n, 3, F) == GR_SUCCESS);
    CHECK(_decfloat_cmp_ui(y, 7000, F) == 0);
    gr_stream_init_str(out);
    CHECK(decfloat_write_sci(out, y, F) == GR_SUCCESS);
    CHECK_STR(out->s, "7e3");
    CHECK(decfloat_via_arb(y, y, _arb_exp, F) == GR_SUCCESS);
    CHECK_STR(decfloat_get_str(y, F), "1.151790051e3040");

    /* operands which are not elements of the context: infinities and nans
       are rejected even where the result would be finite, while any
       precision and exponent is accepted */
    {
        gr_ctx_t Fi;
        decfloat_t inf, nan;
        int c;
        gr_ctx_init_decfloat(Fi, 10, DECIMAL_ALLOW_INF | DECIMAL_ALLOW_NAN);
        decfloat_init(inf, Fi);
        decfloat_init(nan, Fi);
        CHECK(decfloat_pos_inf(inf, Fi) == GR_SUCCESS && decfloat_nan(nan, Fi) == GR_SUCCESS);
        CHECK(decfloat_pos_inf(x, F) == GR_DOMAIN && decfloat_nan(x, F) == GR_UNABLE);
        CHECK(decfloat_inv(x, inf, F) == GR_DOMAIN && decfloat_inv(x, inf, Fi) == GR_SUCCESS && DECFLOAT_IS_ZERO(x));
        CHECK(decfloat_div(x, y, inf, F) == GR_DOMAIN && decfloat_pow_ui(x, inf, 0, F) == GR_DOMAIN && decfloat_pow_ui(x, nan, 0, F) == GR_UNABLE);
        CHECK(decfloat_atan(x, inf, F) == GR_DOMAIN && decfloat_atan(x, inf, Fi) == GR_SUCCESS);
        CHECK(decfloat_cmp(&c, y, inf, F) == GR_DOMAIN && decfloat_cmp(&c, y, inf, Fi) == GR_SUCCESS && c < 0);
        CHECK(decfloat_cmp(&c, y, nan, Fi) == GR_UNABLE && decfloat_sgn(x, inf, F) == GR_DOMAIN);
        CHECK(decfloat_is_integer(inf, F) == T_FALSE && decfloat_is_integer(nan, F) == T_UNKNOWN && decfloat_is_integer(y, F) == T_TRUE);
        CHECK(decfloat_set_si_10exp_si(x, 1, 200, Fi) == GR_SUCCESS);
        decimal_ctx_set_exp_limits(F, -100, 100);
        CHECK(decfloat_set(y, x, F) == GR_UNABLE && decfloat_mul_10exp_si(y, x, -200, F) == GR_SUCCESS && _decfloat_cmp_ui(y, 1, F) == 0);
        CHECK(decfloat_div(y, x, x, F) == GR_SUCCESS && _decfloat_cmp_ui(y, 1, F) == 0);
        CHECK(decfloat_sqrt(y, x, F) == GR_SUCCESS && decfloat_set_round(y, y, 40, DECIMAL_RND_NEAR, F) == GR_SUCCESS);
        decimal_ctx_set_exp_limits(F, WORD_MIN, WORD_MAX);
        decfloat_clear(inf, Fi);
        decfloat_clear(nan, Fi);
        gr_ctx_clear(Fi);
    }

    /* decball */
    CHECK(decball_set_str(X, "[2 +/- 3]", B) == GR_SUCCESS);
    CHECK(_decball_contains_positive(X, B) && _decball_contains_negative(X, B));
    CHECK(_decball_contains_nonpositive(X, B) && _decball_contains_nonnegative(X, B));
    CHECK(decball_set_str(X, "[2 +/- 1]", B) == GR_SUCCESS);
    CHECK(_decball_contains_positive(X, B) && !_decball_contains_negative(X, B));
    CHECK(!_decball_contains_nonpositive(X, B) && _decball_contains_nonnegative(X, B));
    decball_get_rad(u, X, B);
    CHECK(_decmag_cmp_10exp_si(u, 0, B) == 0);
    _decmag_set_10exp_si(u, -3, B);
    CHECK(decball_add_error_decmag(X, u, B) == GR_SUCCESS);
    CHECK(decfloat_set_si_10exp_si(x, -2, -3, F) == GR_SUCCESS);
    CHECK(decball_add_error_decfloat(X, x, B) == GR_SUCCESS);
    CHECK_STR(decball_get_str(X, B), "[2 +/- 1.003]");
    CHECK(decball_via_arb(Y, X, _arb_exp, B) == GR_SUCCESS);
    CHECK(decball_set_str(X, "[7.389056099 +/- 1e-9]", B) == GR_SUCCESS);
    CHECK(_decball_overlaps(X, Y, B));
    CHECK(decball_zero_pm_inf(X, B) == GR_SUCCESS);
    CHECK(!_decball_is_finite(X, B) && _decball_contains_positive(X, B) && _decball_contains_negative(X, B));
    CHECK(gr_fac_ui(X, 20, B) == GR_SUCCESS);
    CHECK_STR(decball_get_str(X, B), "[2432902008000000000 +/- 1.767e8]");

    /* deccfloat */
    CHECK(decfloat_set_si_10exp_si(x, 3, 0, F) == GR_SUCCESS);
    CHECK(decfloat_set_si_10exp_si(y, -4, 0, F) == GR_SUCCESS);
    CHECK(deccfloat_set_decfloat_decfloat(z, x, y, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(z, C), "(3 - 4*I)");
    CHECK(deccfloat_mul_i(w, z, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "(4 + 3*I)");
    CHECK(deccfloat_div_i(w, w, C) == GR_SUCCESS && deccfloat_equal(w, z, C) == T_TRUE);
    CHECK(deccfloat_mul_two(w, z, C) == GR_SUCCESS);
    CHECK(deccfloat_mul_10exp_si(w, w, -1, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "(0.6 - 0.8*I)");
    fmpz_set_si(n, 1);
    CHECK(deccfloat_mul_10exp_fmpz(w, w, n, C) == GR_SUCCESS);
    CHECK(deccfloat_mul_decfloat(w, w, x, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "(18 - 24*I)");
    CHECK(deccfloat_div_decfloat(w, w, y, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "(-4.5 + 6*I)");
    CHECK(deccfloat_set_decfloat(w, y, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "-4");
    CHECK(gr_fac_ui(w, 5, C) == GR_SUCCESS);
    CHECK_STR(deccfloat_get_str(w, C), "120");
    for (i = 0; i < 100; i++)
        GR_IGNORE(deccfloat_randtest_special(w, state, C));

    /* deccball */
    CHECK(decball_set_str(X, "[3 +/- 0.1]", B) == GR_SUCCESS);
    CHECK(decball_set_str(Y, "-4", B) == GR_SUCCESS);
    CHECK(deccball_set_decball_decball(Z, X, Y, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(Z, CB), "([3 +/- 0.1] - 4*I)");
    CHECK(_deccball_contains_deccfloat(Z, z, CB) && !_deccball_contains_deccfloat(Z, w, CB));
    CHECK(deccball_is_integer(Z, CB) == T_FALSE);
    CHECK(deccball_mul_i(W, Z, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(W, CB), "(4 + [3 +/- 0.1]*I)");
    CHECK(deccball_div_i(W, W, CB) == GR_SUCCESS);
    CHECK(deccball_mul_two(W, W, CB) == GR_SUCCESS);
    CHECK(deccball_mul_10exp_si(W, W, 2, CB) == GR_SUCCESS);
    CHECK(deccball_mul_10exp_fmpz(W, W, n, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(W, CB), "([6000 +/- 200] - 8000*I)");
    CHECK(deccball_mul_decball(W, Z, Y, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(W, CB), "([-12 +/- 0.4] + 16*I)");
    CHECK(deccball_div_decball(W, W, Y, CB) == GR_SUCCESS && _deccball_contains(W, Z, CB));
    _decmag_set_10exp_si(u, -2, B);
    CHECK(deccball_add_error_decmag(W, u, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(W, CB), "([3 +/- 0.11] + [-4 +/- 0.01]*I)");
    CHECK(deccball_set_decball(W, Y, CB) == GR_SUCCESS && deccball_is_integer(W, CB) == T_TRUE);
    CHECK(deccball_set_str(W, "([+/- inf] + [3 +/- inf]*I)", CB) == GR_SUCCESS && deccball_is_integer(W, CB) == T_UNKNOWN);
    CHECK(DECMAG_IS_INF(&W->re.rad) && DECMAG_IS_INF(&W->im.rad) && DECFLOAT_IS_ZERO(&W->re.mid) && _decfloat_cmp_ui(&W->im.mid, 3, F) == 0);
    CHECK(deccball_set_str(W, "(inf + 2*I)", CB) != GR_SUCCESS && deccball_set_str(W, "nan", CB) != GR_SUCCESS);
    /* balls are real numbers: no infinite midpoints */
    CHECK(decball_set_str(X, "[+/- inf]", B) == GR_SUCCESS && decball_lower(Y, X, B) == GR_DOMAIN && decball_upper(Y, X, B) == GR_DOMAIN);
    CHECK(decball_abs_upper(Y, X, B) == GR_DOMAIN && decball_abs_lower(Y, X, B) == GR_SUCCESS && decball_is_zero(Y, B) == T_TRUE);
    CHECK(decball_set_str(X, "inf", B) != GR_SUCCESS && decball_set_str(X, "nan", B) != GR_SUCCESS && decball_set_str(X, "1/0", B) != GR_SUCCESS);
    CHECK(decball_set_str(X, "[1 +/- inf]", B) == GR_SUCCESS && decball_mid(Y, X, B) == GR_SUCCESS && decball_is_one(Y, B) == T_TRUE);
    CHECK(decball_rad(Y, X, B) == GR_DOMAIN);
    {
        arb_t a;
        arb_init(a);
        arb_pos_inf(a);
        CHECK(decball_set_arb(X, a, B) == GR_DOMAIN);
        arb_indeterminate(a);
        CHECK(decball_set_arb(X, a, B) == GR_UNABLE);
        arb_zero_pm_inf(a);
        CHECK(decball_set_arb(X, a, B) == GR_SUCCESS && DECMAG_IS_INF(&X->rad));
        arb_clear(a);
    }
    CHECK(gr_fac_ui(W, 5, CB) == GR_SUCCESS && gr_rising_ui(W, W, 2, CB) == GR_SUCCESS);
    CHECK_STR(deccball_get_str(W, CB), "14520");
    fmpq_set_si(q, 1, 2);
    CHECK(gr_gamma_fmpq(W, q, CB) == GR_SUCCESS && gr_gamma_fmpq(X, q, B) == GR_SUCCESS);
    CHECK(_decball_overlaps(X, &W->re, B) && decball_is_zero(&W->im, B) == T_TRUE);

    _decmag_clear(u, B);
    _decmag_clear(v, B);
    decfloat_clear(x, F);
    decfloat_clear(y, F);
    decball_clear(X, B);
    decball_clear(Y, B);
    deccfloat_clear(z, C);
    deccfloat_clear(w, C);
    deccball_clear(Z, CB);
    deccball_clear(W, CB);
    fmpz_clear(n);
    fmpq_clear(q);

    gr_ctx_clear(F);
    gr_ctx_clear(B);
    gr_ctx_clear(C);
    gr_ctx_clear(CB);

    /* digits, limbs and radix integers */
    {
        gr_ctx_t D;
        decfloat_t a, b;
        radix_integer_t m;
        fmpz_t e;
        slong k;

        _gr_ctx_init_decimal(D, DECIMAL_CTX_FLOAT, 3, 30, DECIMAL_RND_NEAR, 0);
        decfloat_init(a, D);
        decfloat_init(b, D);
        radix_integer_init(m, DECIMAL_CTX_RADIX(D));
        fmpz_init(e);

        CHECK(decfloat_set_str(a, "-12345.678", D) == GR_SUCCESS);
        CHECK(decfloat_digits(a, D) == 8 && decfloat_limbs(a, D) == 3);
        CHECK(decfloat_get_digit_si(a, 4, D) == 1 && decfloat_get_digit_si(a, 0, D) == 5 && decfloat_get_digit_si(a, -3, D) == 8);
        CHECK(decfloat_get_digit_si(a, 5, D) == 0 && decfloat_get_digit_si(a, -4, D) == 0 && decfloat_get_digit_si(a, -1000, D) == 0);
        CHECK(decfloat_get_limb_si(a, 1, D) == 12 && decfloat_get_limb_si(a, 0, D) == 345 && decfloat_get_limb_si(a, -1, D) == 678);
        CHECK(decfloat_set_digit_si(b, a, 7, 4, D) == GR_SUCCESS);
        CHECK_STR(decfloat_get_str(b, D), "-40012345.678");
        CHECK(decfloat_set_digit_si(b, b, -9, 1, D) == GR_SUCCESS);
        CHECK_STR(decfloat_get_str(b, D), "-40012345.678000001");
        CHECK(decfloat_set_limb_si(b, a, -1, 0, D) == GR_SUCCESS);
        CHECK_STR(decfloat_get_str(b, D), "-12345");
        CHECK(decfloat_set_limb_si(b, b, 1, 0, D) == GR_SUCCESS && decfloat_set_limb_si(b, b, 0, 0, D) == GR_SUCCESS && DECFLOAT_IS_ZERO(b));
        CHECK(decfloat_set_limb_si(b, b, 2, 7, D) == GR_SUCCESS);
        CHECK_STR(decfloat_get_str(b, D), "7000000");
        CHECK(decfloat_set_limb_si(b, a, 0, 1000, D) == GR_DOMAIN && decfloat_set_digit_si(b, a, 0, 10, D) == GR_DOMAIN);

        CHECK(decfloat_get_radix_integer(m, a, D) == GR_DOMAIN);
        CHECK(decfloat_get_radix_integer_Bexp_fmpz(m, e, a, D) == GR_SUCCESS && fmpz_equal_si(e, -1));
        CHECK(decfloat_set_radix_integer_Bexp_fmpz(b, m, e, D) == GR_SUCCESS && decfloat_equal(a, b, D) == T_TRUE);
        CHECK(decfloat_set_str(a, "-12e6", D) == GR_SUCCESS && decfloat_get_radix_integer(m, a, D) == GR_SUCCESS);
        CHECK(radix_integer_size_limbs(m, DECIMAL_CTX_RADIX(D)) == 3 && radix_integer_get_limb(m, 2, DECIMAL_CTX_RADIX(D)) == 12 && m->size < 0);
        CHECK(decfloat_set_radix_integer(b, m, D) == GR_SUCCESS && decfloat_equal(a, b, D) == T_TRUE);
        decimal_ctx_set_prec(D, 1);
        CHECK(decfloat_set_radix_integer(b, m, D) == GR_SUCCESS);
        CHECK_STR(decfloat_get_str(b, D), "-10000000");

        /* random round trips through the digits */
        decimal_ctx_set_prec(D, 30);
        for (i = 0; i < 100; i++)
        {
            GR_MUST_SUCCEED(decfloat_randtest(a, state, D));
            if (!DECFLOAT_IS_FINITE(a) || COEFF_IS_MPZ(a->exp) || FLINT_ABS(a->exp) > 100)
                continue;
            GR_MUST_SUCCEED(decfloat_zero(b, D));
            for (k = -400; k <= 400; k++)
                if (decfloat_get_digit_si(a, k, D) != 0)
                    GR_MUST_SUCCEED(decfloat_set_digit_si(b, b, k, decfloat_get_digit_si(a, k, D), D));
            if (_decfloat_sgn(a, D) < 0)
                GR_MUST_SUCCEED(decfloat_neg(b, b, D));
            CHECK(decfloat_equal(a, b, D) == T_TRUE);
        }

        decfloat_clear(a, D);
        decfloat_clear(b, D);
        radix_integer_clear(m, DECIMAL_CTX_RADIX(D));
        fmpz_clear(e);
        gr_ctx_clear(D);
    }

    /* long operands: the arb enclosure only uses the leading limbs, and
       exact values are detected from the digits */
    {
        gr_ctx_t X, R;
        decfloat_t a, b;
        arb_t t;
        fmpq_t qa;
        char * str;
        slong L = 3000, prec_bits;

        gr_ctx_init_decfloat(X, DECIMAL_PREC_EXACT, 0);
        gr_ctx_init_decfloat(R, 20, 0);
        decimal_ctx_set_rnd(R, DECIMAL_RND_DOWN);
        decfloat_init(a, X);
        decfloat_init(b, X);
        arb_init(t);
        fmpq_init(qa);
        str = flint_malloc(L + 20);

        for (i = 0; i < 20; i++)
        {
            slong j;
            str[0] = '1' + n_randint(state, 9);
            str[1] = '.';
            for (j = 2; j < L; j++)
                str[j] = (n_randint(state, 4) == 0) ? '0' + n_randint(state, 10) : '0';
            str[L] = '1';
            flint_sprintf(str + L + 1, "e%wd", (slong) n_randint(state, 2000) - 1000);
            GR_MUST_SUCCEED(decfloat_set_str(a, str, X));
            if (n_randint(state, 2))
                GR_MUST_SUCCEED(decfloat_neg(a, a, X));
            prec_bits = 2 + n_randint(state, 300);
            GR_MUST_SUCCEED(decfloat_get_arb(t, a, prec_bits, X));
            GR_MUST_SUCCEED(decfloat_get_fmpq(qa, a, X));
            CHECK(arb_contains_fmpq(t, qa) && arb_rel_accuracy_bits(t) >= prec_bits - 5);
        }

        /* 1 + 10^-L: elementary functions go up to the length of the operand */
        str[0] = '1'; str[1] = '.';
        memset(str + 2, '0', L);
        strcpy(str + 2 + L, "1");
        GR_MUST_SUCCEED(decfloat_set_str(a, str, X));
        CHECK(gr_log2(b, a, R) == GR_SUCCESS);
        {
            char * t = decfloat_get_str(b, R);
            CHECK(t[0] == '1' && _decfloat_sci_exp_clamped(b, R) == -L - 1);
            flint_free(t);
        }

        /* 2^-1000 has 699 digits: log2 detects it via digit tests and a modular check */
        GR_MUST_SUCCEED(decfloat_set_ui(a, 2, X));
        GR_MUST_SUCCEED(decfloat_pow_si(a, a, -1000, X));
        CHECK(gr_log2(b, a, R) == GR_SUCCESS && _decfloat_cmp_si(b, -1000, R) == 0);
        GR_MUST_SUCCEED(decfloat_set_str(b, "1e-1000", X));
        GR_MUST_SUCCEED(decfloat_add(a, a, b, X));
        CHECK(gr_log2(b, a, R) == GR_SUCCESS && _decfloat_cmp_si(b, -1000, R) > 0);

        /* sin(pi x) for a long half-integer, zeta at a long even negative integer */
        memset(str, '7', L);
        strcpy(str + L, ".5");
        GR_MUST_SUCCEED(decfloat_set_str(a, str, X));
        CHECK(gr_sin_pi(b, a, R) == GR_SUCCESS && _decfloat_cmp_si(b, -1, R) == 0);
        str[L - 1] = '4';
        str[L] = '\0';
        GR_MUST_SUCCEED(decfloat_set_str(a, str, X));
        GR_MUST_SUCCEED(decfloat_neg(a, a, X));
        CHECK(gr_zeta(b, a, R) == GR_SUCCESS && DECFLOAT_IS_ZERO(b));

        decfloat_clear(a, X);
        decfloat_clear(b, X);
        arb_clear(t);
        fmpq_clear(qa);
        flint_free(str);
        gr_ctx_clear(X);
        gr_ctx_clear(R);
    }

    TEST_FUNCTION_END(state);
}

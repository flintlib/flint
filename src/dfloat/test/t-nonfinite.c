/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include <string.h>
#include "test_helpers.h"
#include "double_extras.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "dfloat.h"

/* The conventions for nonfinite balls and domains (see the section
   "Nonfinite values and domains" of the documentation), through the
   generic rings: every form of the whole line W, the status codes of
   the functions on and off their domains, the normal form of W in
   results, and the documented examples. */

#define NF_CHECK(cond, what) \
    do { if (!(cond)) { flint_printf("FAIL: %s: %s (line %d, d%db flags %d, form %d)\n", \
        what, #cond, __LINE__, n, flags, form); flint_abort(); } } while (0)

#define NF_ELEM 5

/* ball element from components (n doubles) and radius */
static void
_nf_set(double * x, int n, double d0, double d1, double rad)
{
    int i;
    for (i = 0; i < n; i++)
        x[i] = 0.0;
    x[0] = d0;
    if (n > 1)
        x[1] = d1;
    x[n] = rad;
}

/* the whole line in normal form */
static int
_nf_is_normal_W(const double * x, int n)
{
    int i;
    for (i = 0; i < n; i++)
        if (x[i] != 0.0)
            return 0;
    return x[n] == D_INF;
}

/* finite, [0 +/- b] with b >= c (c > 0) or the whole line */
static int
_nf_covers_zero_pm(const double * x, int n, double c)
{
    int i;
    if (_nf_is_normal_W(x, n))
        return 1;
    for (i = 0; i < n; i++)
        if (x[i] != 0.0)
            return 0;
    return x[n] >= c && x[n] <= 2.0 * c;
}

/* the forms of the whole line */
#define NF_FORMS 10
static void
_nf_W(double * x, int n, int form)
{
    switch (form)
    {
        case 0: _nf_set(x, n, 0.0, 0.0, D_INF); break;
        case 1: _nf_set(x, n, 5.0, 0.0, D_INF); break;
        case 2: _nf_set(x, n, D_INF, 0.0, 0.0); break;
        case 3: _nf_set(x, n, -D_INF, 0.0, 1.0); break;
        case 4: _nf_set(x, n, D_NAN, 0.0, 0.0); break;
        case 5: _nf_set(x, n, 1.0, 0.0, D_NAN); break;
        case 6: _nf_set(x, n, D_NAN, 0.0, D_INF); break;
        case 7: _nf_set(x, n, D_INF, 0.0, D_INF); break;
        case 8: _nf_set(x, n, 0.0, 0.0, D_NAN); break;
        default:
            /* a nonfinite tail (n >= 2), else a nonfinite head */
            if (n >= 2) _nf_set(x, n, 1.0, D_NAN, 0.0);
            else _nf_set(x, n, D_NAN, 0.0, 0.0);
            break;
    }
}

TEST_FUNCTION_START(nonfinite, state)
{
    static const int flagsets[3] = { DFLOAT_BALL, DFLOAT_BALL | DFLOAT_STRONG, DFLOAT_BALL | DFLOAT_FAST };
    int n, fl, flags, form;
    gr_ctx_t ctx;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (n = 1; n <= 4; n++)
    for (fl = 0; fl < 3; fl++)
    {
        double x[NF_ELEM], y[NF_ELEM], z[NF_ELEM], w[NF_ELEM], one[NF_ELEM], zero[NF_ELEM], three[NF_ELEM];
        char * s0, * s1;
        int cmp;

        flags = flagsets[fl];
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, flags));
        form = -1;

        _nf_set(one, n, 1.0, 0.0, 0.0);
        _nf_set(zero, n, 0.0, 0.0, 0.0);
        _nf_set(three, n, 3.0, 0.0, 0.0);

        /* the documented examples */
        _nf_W(x, n, 0);
        NF_CHECK(gr_div(z, x, three, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "W / 3");
        NF_CHECK(gr_div(z, one, zero, ctx) == GR_DOMAIN, "1 / 0");
        NF_CHECK(gr_div(z, x, zero, ctx) == GR_DOMAIN, "W / 0");
        _nf_set(y, n, 3.0, 0.0, 2.0);
        NF_CHECK(gr_div(z, x, y, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "W / [3 +/- 2]");
        _nf_set(y, n, 1.0, 0.0, 2.0);
        NF_CHECK(gr_div(z, one, y, ctx) == GR_UNABLE, "1 / [1 +/- 2]");
        NF_CHECK(gr_div(z, one, x, ctx) == GR_UNABLE, "1 / W");
        NF_CHECK(gr_exp(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "exp(W)");
        _nf_set(y, n, D_INF, 0.0, 0.0);
        NF_CHECK(gr_sin(z, y, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.0), "sin([inf +/- 0])");
        _nf_set(y, n, -1.0, 0.0, 0.0);
        NF_CHECK(gr_sqrt(z, y, ctx) == GR_DOMAIN, "sqrt(-1)");
        _nf_set(y, n, 0.0, 0.0, 1.0);
        NF_CHECK(gr_sqrt(z, y, ctx) == GR_UNABLE, "sqrt([0 +/- 1])");
        NF_CHECK(gr_log(z, zero, ctx) == GR_DOMAIN, "log(0)");
        _nf_set(y, n, -2.0, 0.0, 1.0);
        NF_CHECK(gr_log(z, y, ctx) == GR_DOMAIN, "log([-2 +/- 1])");
        NF_CHECK(gr_log(z, x, ctx) == GR_UNABLE, "log(W)");
        NF_CHECK(gr_asin(z, one, ctx) == GR_SUCCESS && fabs(z[0] - 1.5707963267948966) < 1e-15, "asin(1)");
        _nf_set(y, n, 2.0, 0.0, 0.5);
        NF_CHECK(gr_asin(z, y, ctx) == GR_DOMAIN, "asin([2 +/- 0.5])");
        _nf_set(y, n, -0.5, 0.0, 0.0);
        NF_CHECK(gr_pow(z, zero, y, ctx) == GR_DOMAIN, "0^(-1/2)");
        _nf_set(w, n, -2.0, 0.0, 0.0);
        _nf_set(y, n, 0.5, 0.0, 0.0);
        NF_CHECK(gr_pow(z, w, y, ctx) == GR_DOMAIN, "(-2)^(1/2)");
        NF_CHECK(gr_pow(z, w, x, ctx) == GR_UNABLE, "(-2)^W");
        NF_CHECK(gr_pow(z, x, three, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "W^3");
        NF_CHECK(gr_pow(z, x, zero, ctx) == GR_SUCCESS && gr_is_one(z, ctx) == T_TRUE, "W^0");
        NF_CHECK(gr_mul(z, x, zero, ctx) == GR_SUCCESS && gr_is_zero(z, ctx) == T_TRUE, "W * 0");

        gr_get_str(&s0, x, ctx);

        for (form = 0; form < NF_FORMS; form++)
        {
            double sn[NF_ELEM], cs[NF_ELEM];
            _nf_W(x, n, form);

            /* entire functions: success, W or the bounded range */
            NF_CHECK(gr_add(z, x, one, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "add");
            NF_CHECK(gr_sub(z, one, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "sub");
            NF_CHECK(gr_mul(z, x, three, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "mul");
            NF_CHECK(gr_mul(z, x, zero, ctx) == GR_SUCCESS && gr_is_zero(z, ctx) == T_TRUE, "mul by 0");
            NF_CHECK(gr_sqr(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "sqr");
            NF_CHECK(gr_neg(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "neg");
            NF_CHECK(gr_abs(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "abs");
            NF_CHECK(gr_set(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "set");
            NF_CHECK(gr_mul_2exp_si(z, x, -3, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "mul_2exp_si");
            NF_CHECK(gr_exp(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "exp");
            NF_CHECK(gr_expm1(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "expm1");
            NF_CHECK(gr_sinh(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "sinh");
            NF_CHECK(gr_cosh(z, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "cosh");
            NF_CHECK(gr_sin(z, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.0), "sin");
            NF_CHECK(gr_cos(z, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.0), "cos");
            NF_CHECK(gr_sin_cos(sn, cs, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(sn, n, 1.0)
                && _nf_covers_zero_pm(cs, n, 1.0), "sin_cos");
            NF_CHECK(gr_tanh(z, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.0), "tanh");
            NF_CHECK(gr_atan(z, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.5707963267948966), "atan");
            NF_CHECK(gr_atan2(z, x, one, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 3.141592653589793), "atan2");
            NF_CHECK(gr_sgn(z, x, ctx) == GR_SUCCESS && _nf_covers_zero_pm(z, n, 1.0), "sgn");
            NF_CHECK(gr_div(z, x, three, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "div (dividend)");

            /* not entire: GR_UNABLE (W contains points on both sides) */
            NF_CHECK(gr_div(z, one, x, ctx) == GR_UNABLE, "div");
            NF_CHECK(gr_div(z, x, x, ctx) == GR_UNABLE, "div aliased");
            NF_CHECK(gr_div(z, x, zero, ctx) == GR_DOMAIN, "div by 0");
            NF_CHECK(gr_inv(z, x, ctx) == GR_UNABLE, "inv");
            NF_CHECK(gr_sqrt(z, x, ctx) == GR_UNABLE, "sqrt");
            NF_CHECK(gr_rsqrt(z, x, ctx) == GR_UNABLE, "rsqrt");
            NF_CHECK(gr_log(z, x, ctx) == GR_UNABLE, "log");
            NF_CHECK(gr_log1p(z, x, ctx) == GR_UNABLE, "log1p");
            NF_CHECK(gr_tan(z, x, ctx) == GR_UNABLE, "tan");
            NF_CHECK(gr_asin(z, x, ctx) == GR_UNABLE, "asin");
            NF_CHECK(gr_acos(z, x, ctx) == GR_UNABLE, "acos");
            _nf_set(y, n, 0.5, 0.0, 0.0);
            NF_CHECK(gr_pow(z, x, y, ctx) == GR_UNABLE, "pow (base)");
            NF_CHECK(gr_pow(z, three, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "pow (exponent, base > 0)");
            NF_CHECK(gr_pow(z, one, x, ctx) == GR_SUCCESS && gr_is_one(z, ctx) == T_TRUE, "1^W");
            NF_CHECK(gr_pow(z, zero, x, ctx) == GR_UNABLE, "0^W");
            NF_CHECK(gr_pow_si(z, x, -2, ctx) == GR_UNABLE, "pow_si (negative)");
            NF_CHECK(gr_pow_si(z, x, 2, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "pow_si (positive)");

            /* predicates and comparisons */
            NF_CHECK(gr_is_zero(x, ctx) == T_UNKNOWN, "is_zero");
            NF_CHECK(gr_is_one(x, ctx) == T_UNKNOWN, "is_one");
            NF_CHECK(gr_is_neg_one(x, ctx) == T_UNKNOWN, "is_neg_one");
            NF_CHECK(gr_is_invertible(x, ctx) == T_UNKNOWN, "is_invertible");
            NF_CHECK(gr_equal(x, x, ctx) == T_UNKNOWN, "equal");
            NF_CHECK(gr_equal(x, one, ctx) == T_UNKNOWN, "equal");
            NF_CHECK(gr_cmp(&cmp, x, one, ctx) == GR_UNABLE, "cmp");
            NF_CHECK(gr_cmpabs(&cmp, one, x, ctx) == GR_UNABLE, "cmpabs");

            /* conversions */
            {
                fmpz_t e;
                slong si;
                double d;
                arb_t a;
                double m[NF_ELEM], r[NF_ELEM];
                gr_ctx_t RR;
                fmpz_init(e);
                arb_init(a);
                gr_ctx_init_real_arb(RR, 64);
                NF_CHECK(gr_get_fmpz(e, x, ctx) == GR_UNABLE, "get_fmpz");
                NF_CHECK(gr_get_si(&si, x, ctx) == GR_UNABLE, "get_si");
                NF_CHECK(gr_get_d(&d, x, ctx) == GR_SUCCESS && d == 0.0, "get_d");
                NF_CHECK(gr_set_other(a, x, ctx, RR) == GR_SUCCESS && !arb_is_finite(a)
                    && arb_contains_si(a, 0) && arb_contains_si(a, -123456), "to arb");
                NF_CHECK(gr_get_interval_mid_rad(m, r, x, ctx) == GR_UNABLE, "get_interval_mid_rad");
                NF_CHECK(gr_set_interval_mid_rad(z, x, one, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "set_interval_mid_rad");
                NF_CHECK(gr_set_interval_mid_rad(z, one, x, ctx) == GR_SUCCESS && _nf_is_normal_W(z, n), "set_interval_mid_rad");
                gr_get_str(&s1, x, ctx);
                NF_CHECK(strcmp(s0, s1) == 0, "printing");
                flint_free(s1);
                gr_ctx_clear(RR);
                arb_clear(a);
                fmpz_clear(e);
            }
            NF_CHECK(gr_set_d(z, x[0], ctx) == (isfinite(x[0]) ? GR_SUCCESS : GR_DOMAIN), "set_d");

            /* vectors, long enough for the SIMD blocks, with W in one
               lane: the other lanes are unaffected */
            {
                double X[12 * NF_ELEM], Y[12 * NF_ELEM], Z[12 * NF_ELEM], Z2[12 * NF_ELEM], R[NF_ELEM];
                int i, len = 12, k = 5, st;
                for (i = 0; i < len; i++)
                {
                    _nf_set(X + i * (n + 1), n, 0.25 + i, 0.0, 0.0);
                    _nf_set(Y + i * (n + 1), n, 2.0 - 0.125 * i, 0.0, 0.0);
                }
                memcpy(X + k * (n + 1), x, sizeof(double) * (n + 1));

#define NF_VEC_OK(V, what) \
                for (i = 0; i < len; i++) \
                    NF_CHECK((i == k) == _nf_is_normal_W(V + i * (n + 1), n), what);
                NF_CHECK(_gr_vec_add(Z, X, Y, len, ctx) == GR_SUCCESS, "vec_add"); NF_VEC_OK(Z, "vec_add")
                NF_CHECK(_gr_vec_sub(Z, Y, X, len, ctx) == GR_SUCCESS, "vec_sub"); NF_VEC_OK(Z, "vec_sub")
                NF_CHECK(_gr_vec_mul(Z, X, Y, len, ctx) == GR_SUCCESS, "vec_mul"); NF_VEC_OK(Z, "vec_mul")
                NF_CHECK(_gr_vec_div(Z, X, Y, len, ctx) == GR_SUCCESS, "vec_div"); NF_VEC_OK(Z, "vec_div")
                NF_CHECK(_gr_vec_mul_scalar(Z, X, len, three, ctx) == GR_SUCCESS, "vec_mul_scalar"); NF_VEC_OK(Z, "vec_mul_scalar")
                memcpy(Z, Y, sizeof(Y));
                NF_CHECK(_gr_vec_addmul_scalar(Z, X, len, three, ctx) == GR_SUCCESS, "vec_addmul_scalar"); NF_VEC_OK(Z, "vec_addmul_scalar")
                NF_CHECK(_gr_vec_exp(Z, X, len, ctx) == GR_SUCCESS, "vec_exp"); NF_VEC_OK(Z, "vec_exp")
                NF_CHECK(_gr_vec_sin_cos(Z, Z2, X, len, ctx) == GR_SUCCESS, "vec_sin_cos");
                NF_CHECK(_nf_covers_zero_pm(Z + k * (n + 1), n, 1.0) && _nf_covers_zero_pm(Z2 + k * (n + 1), n, 1.0), "vec_sin_cos");
                st = _gr_vec_div(Z, Y, X, len, ctx);
                NF_CHECK(st == GR_UNABLE, "vec_div (divisor)");
                NF_CHECK(_gr_vec_log(Z, X, len, ctx) == GR_UNABLE, "vec_log");
                NF_CHECK(_gr_vec_sqrt(Z, X, len, ctx) == GR_UNABLE, "vec_sqrt");
                NF_CHECK(_gr_vec_rsqrt(Z, X, len, ctx) == GR_UNABLE, "vec_rsqrt");

                /* dot products: W times an exact zero is zero */
                NF_CHECK(_gr_vec_dot(R, NULL, 0, X, Y, len, ctx) == GR_SUCCESS && _nf_is_normal_W(R, n), "vec_dot");
                NF_CHECK(_gr_vec_dot_rev(R, NULL, 1, X, Y, len, ctx) == GR_SUCCESS && _nf_is_normal_W(R, n), "vec_dot_rev");
                _nf_set(Y + k * (n + 1), n, 0.0, 0.0, 0.0);
                NF_CHECK(_gr_vec_dot(R, NULL, 0, X, Y, len, ctx) == GR_SUCCESS && !_nf_is_normal_W(R, n) && R[n] < 1e-10, "vec_dot (W * 0)");
                NF_CHECK(_gr_vec_dot(R, x, 0, Y, Y, 3, ctx) == GR_SUCCESS && _nf_is_normal_W(R, n), "vec_dot (initial W)");
                NF_CHECK(_gr_vec_dot(R, x, 0, Y, Y, 0, ctx) == GR_SUCCESS && _nf_is_normal_W(R, n), "vec_dot (initial W, len 0)");
#undef NF_VEC_OK
            }
        }

        flint_free(s0);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}

#undef NF_CHECK
#undef NF_ELEM
#undef NF_FORMS

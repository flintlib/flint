/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "decimal.h"
#include "mag.h"
#include "fmpz_extras.h"
#include "gr.h"

/* ------------------------------------------------------------------------- */
/*    Special values                                                         */
/* ------------------------------------------------------------------------- */

void _decfloat_pos_inf(decfloat_t res)
{
    res->m.size = 0;
    fmpz_set_si(&res->exp, DECFLOAT_EXP_POS_INF);
}

void _decfloat_neg_inf(decfloat_t res)
{
    res->m.size = 0;
    fmpz_set_si(&res->exp, DECFLOAT_EXP_NEG_INF);
}

void _decfloat_nan(decfloat_t res)
{
    res->m.size = 0;
    fmpz_set_si(&res->exp, DECFLOAT_EXP_NAN);
}

int
decfloat_pos_inf(decfloat_t res, gr_ctx_t ctx)
{
    if (!DECIMAL_CTX_ALLOW_INF(ctx))
        return GR_DOMAIN;
    _decfloat_pos_inf(res);
    return GR_SUCCESS;
}

int
decfloat_neg_inf(decfloat_t res, gr_ctx_t ctx)
{
    if (!DECIMAL_CTX_ALLOW_INF(ctx))
        return GR_DOMAIN;
    _decfloat_neg_inf(res);
    return GR_SUCCESS;
}

int
decfloat_nan(decfloat_t res, gr_ctx_t ctx)
{
    if (!DECIMAL_CTX_ALLOW_NAN(ctx))
        return GR_UNABLE;
    _decfloat_nan(res);
    return GR_SUCCESS;
}

int
decfloat_one(decfloat_t res, gr_ctx_t ctx)
{
    radix_integer_fit_limbs(&res->m, 1, DECIMAL_CTX_RADIX(ctx))[0] = 1;
    res->m.size = 1;
    fmpz_zero(&res->exp);
    return _decfloat_finalize(res, ctx);
}

int
decfloat_neg_one(decfloat_t res, gr_ctx_t ctx)
{
    radix_integer_fit_limbs(&res->m, 1, DECIMAL_CTX_RADIX(ctx))[0] = 1;
    res->m.size = -1;
    fmpz_zero(&res->exp);
    return _decfloat_finalize(res, ctx);
}

/* ------------------------------------------------------------------------- */
/*    Predicates                                                             */
/* ------------------------------------------------------------------------- */

truth_t
decfloat_is_zero(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x))
        return T_UNKNOWN;
    return DECFLOAT_IS_ZERO(x) ? T_TRUE : T_FALSE;
}

truth_t
decfloat_is_one(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x))
        return T_UNKNOWN;
    return (x->m.size == 1 && x->m.d[0] == 1 && fmpz_is_zero(&x->exp)) ? T_TRUE : T_FALSE;
}

truth_t
decfloat_is_neg_one(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x))
        return T_UNKNOWN;
    return (x->m.size == -1 && x->m.d[0] == 1 && fmpz_is_zero(&x->exp)) ? T_TRUE : T_FALSE;
}

truth_t
decfloat_is_integer(const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x))
        return T_UNKNOWN;
    return _decfloat_is_int(x, ctx) ? T_TRUE : T_FALSE;
}

truth_t
decfloat_equal(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return T_UNKNOWN;

    if (x->m.size != y->m.size)
        return T_FALSE;

    if (!fmpz_equal(&x->exp, &y->exp))
        return T_FALSE;

    if (x->m.size == 0)
        return T_TRUE;

    return flint_mpn_equal_p(x->m.d, y->m.d, FLINT_ABS(x->m.size)) ? T_TRUE : T_FALSE;
}

int
_decfloat_sgn(const decfloat_t x, gr_ctx_t ctx)
{
    if (x->m.size != 0)
        return (x->m.size > 0) ? 1 : -1;
    if (DECFLOAT_IS_POS_INF(x))
        return 1;
    if (DECFLOAT_IS_NEG_INF(x))
        return -1;
    return 0;
}

/* Compare |x| and |y| for nonzero finite x, y. */
int
_decfloat_cmpabs_finite(const decfloat_t x, const decfloat_t y)
{
    slong xn = FLINT_ABS(x->m.size);
    slong yn = FLINT_ABS(y->m.size);
    slong shift, i;

    /* compare top positions: (vx + xn) vs (vy + yn) */
    shift = _fmpz_sub_small(&x->exp, &y->exp);   /* saturating */

    if (shift > WORD_MAX / 4)
        return 1;
    if (shift < -WORD_MAX / 4)
        return -1;

    /* top_x - top_y = shift + xn - yn */
    shift += xn - yn;

    if (shift != 0)
        return (shift > 0) ? 1 : -1;

    /* aligned tops: compare limbwise from the top, missing limbs are zero */
    for (i = 0; i < FLINT_MAX(xn, yn); i++)
    {
        ulong a = (i < xn) ? x->m.d[xn - 1 - i] : 0;
        ulong b = (i < yn) ? y->m.d[yn - 1 - i] : 0;

        if (a != b)
            return (a < b) ? -1 : 1;
    }

    return 0;
}

int
_decfloat_cmp(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    int sx, sy;

    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y))
    {
        if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
            return 0;

        sx = _decfloat_sgn(x, ctx);
        sy = _decfloat_sgn(y, ctx);

        if (DECFLOAT_IS_INF(x) || DECFLOAT_IS_INF(y))
        {
            if (DECFLOAT_IS_INF(x) && DECFLOAT_IS_INF(y))
                return (sx < sy) ? -1 : (sx > sy);
            if (DECFLOAT_IS_INF(x))
                return sx;
            return -sy;
        }

        /* one of them is zero */
        return (sx < sy) ? -1 : (sx > sy);
    }

    sx = (x->m.size > 0) ? 1 : -1;
    sy = (y->m.size > 0) ? 1 : -1;

    if (sx != sy)
        return sx;

    return sx * _decfloat_cmpabs_finite(x, y);
}

int
_decfloat_cmpabs(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y))
    {
        int ix, iy;

        if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
            return 0;

        ix = DECFLOAT_IS_INF(x) ? 2 : (DECFLOAT_IS_ZERO(x) ? 0 : 1);
        iy = DECFLOAT_IS_INF(y) ? 2 : (DECFLOAT_IS_ZERO(y) ? 0 : 1);

        return (ix < iy) ? -1 : (ix > iy);
    }

    return _decfloat_cmpabs_finite(x, y);
}

int
_decfloat_cmp_ui(const decfloat_t x, ulong c, gr_ctx_t ctx)
{
    decfloat_t t;
    int r;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_ui(t, c, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
    r = _decfloat_cmp(x, t, ctx);
    decfloat_clear(t, ctx);
    return r;
}

int
_decfloat_cmp_si(const decfloat_t x, slong c, gr_ctx_t ctx)
{
    decfloat_t t;
    int r;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_si(t, c, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
    r = _decfloat_cmp(x, t, ctx);
    decfloat_clear(t, ctx);
    return r;
}

int
_decfloat_cmpabs_ui(const decfloat_t x, ulong c, gr_ctx_t ctx)
{
    decfloat_t t;
    int r;
    decfloat_init(t, ctx);
    GR_MUST_SUCCEED(decfloat_set_round_ui(t, c, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));
    r = _decfloat_cmpabs(x, t, ctx);
    decfloat_clear(t, ctx);
    return r;
}

int
decfloat_cmp(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    *res = 0;
    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return GR_UNABLE;
    DECFLOAT_CHECK_OPERAND(x, ctx);
    DECFLOAT_CHECK_OPERAND(y, ctx);
    *res = _decfloat_cmp(x, y, ctx);
    return GR_SUCCESS;
}

int
decfloat_cmpabs(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    *res = 0;
    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return GR_UNABLE;
    DECFLOAT_CHECK_OPERAND(x, ctx);
    DECFLOAT_CHECK_OPERAND(y, ctx);
    *res = _decfloat_cmpabs(x, y, ctx);
    return GR_SUCCESS;
}

int
decfloat_sgn(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x))
        return decfloat_nan(res, ctx);
    DECFLOAT_CHECK_OPERAND(x, ctx);
    return decfloat_set_si(res, _decfloat_sgn(x, ctx), ctx);
}

/* ------------------------------------------------------------------------- */
/*    Sign manipulation                                                      */
/* ------------------------------------------------------------------------- */

int
decfloat_neg_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    int status;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        /* internal exact negations (NOLIMITS) also copy special values
           that the context does not admit */
        if (rnd & DECIMAL_RND_NOLIMITS)
        {
            if (DECFLOAT_IS_POS_INF(x))
                _decfloat_neg_inf(res);
            else if (DECFLOAT_IS_NEG_INF(x))
                _decfloat_pos_inf(res);
            else if (DECFLOAT_IS_NAN(x))
                _decfloat_nan(res);
            else
                decfloat_zero(res, ctx);
            return GR_SUCCESS;
        }
        if (DECFLOAT_IS_POS_INF(x))
            return decfloat_neg_inf(res, ctx);
        if (DECFLOAT_IS_NEG_INF(x))
            return decfloat_pos_inf(res, ctx);
        return decfloat_set_round(res, x, prec, rnd, ctx);
    }

    /* negating flips floor/ceil */
    rnd = DECIMAL_RND_NEGATE(rnd);

    status = decfloat_set_round(res, x, prec, rnd, ctx);
    res->m.size = -res->m.size;
    return status;
}

int
decfloat_neg(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_neg_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_abs_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    if (x->m.size < 0)
        return decfloat_neg_round(res, x, prec, rnd, ctx);
    if (DECFLOAT_IS_NEG_INF(x))
        return decfloat_neg_round(res, x, prec, rnd, ctx);
    return decfloat_set_round(res, x, prec, rnd, ctx);
}

int
decfloat_abs(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_abs_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

/* The remaining functions are not performance-critical. */
PUSH_OPTIONS
OPTIMIZE_OSIZE

/* ------------------------------------------------------------------------- */
/*    Size information                                                       */
/* ------------------------------------------------------------------------- */

slong
_decfloat_mant_digits(const decfloat_t x, gr_ctx_t ctx)
{
    slong n = FLINT_ABS(x->m.size);
    if (n == 0)
        return 0;
    return (n - 1) * DECIMAL_CTX_E(ctx) + _radix_size_digits_1(x->m.d[n - 1], DECIMAL_CTX_RADIX(ctx));
}

slong
decfloat_digits(const decfloat_t x, gr_ctx_t ctx)
{
    slong n = FLINT_ABS(x->m.size);
    if (n == 0)
        return 0;
    return _decfloat_mant_digits(x, ctx) - _radix_valuation_digits_1(x->m.d[0], DECIMAL_CTX_RADIX(ctx));
}

slong
decfloat_limbs(const decfloat_t x, gr_ctx_t ctx)
{
    return FLINT_ABS(x->m.size);
}

/* v with |x| = M 10^v, 10 not dividing M, clamped to +/- DECFLOAT_EXP_CLAMP;
   x finite and nonzero */
slong
_decfloat_val10_clamped(const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong v = _radix_valuation_digits_1(x->m.d[0], DECIMAL_CTX_RADIX(ctx));

    FLINT_ASSERT(x->m.size != 0);

    if (COEFF_IS_MPZ(x->exp) || x->exp > DECFLOAT_EXP_CLAMP / e || x->exp < -DECFLOAT_EXP_CLAMP / e)
        return (fmpz_sgn(&x->exp) > 0) ? DECFLOAT_EXP_CLAMP : -DECFLOAT_EXP_CLAMP;

    return x->exp * e + v;
}

/* index of the limb of x at B^k, or -1 if out of range */
static slong
_decfloat_limb_index(const decfloat_t x, slong k)
{
    slong j;

    if (COEFF_IS_MPZ(x->exp))
        return -1;
    /* |k| and |exp| are both below 2^(FLINT_BITS - 2) */
    j = k - x->exp;
    return (j >= 0 && j < FLINT_ABS(x->m.size)) ? j : -1;
}

ulong
decfloat_get_limb_si(const decfloat_t x, slong k, gr_ctx_t ctx)
{
    slong j;

    if (DECFLOAT_IS_SPECIAL(x))
        return 0;

    j = _decfloat_limb_index(x, k);
    return (j >= 0) ? x->m.d[j] : 0;
}

ulong
decfloat_get_digit_si(const decfloat_t x, slong k, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong q = (k >= 0) ? k / e : -((-k + e - 1) / e);    /* floor(k / e) */
    ulong c = decfloat_get_limb_si(x, q, ctx);

    return (c / n_pow(10, k - q * e)) % 10;
}

int
decfloat_set_limb_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong n, lo, hi, span, j, a, b;
    int negative;
    nn_ptr t, d;
    TMP_INIT;

    if (c >= DECIMAL_CTX_B(ctx) || (DECFLOAT_IS_SPECIAL(x) && !DECFLOAT_IS_ZERO(x)))
        return GR_DOMAIN;

    if (DECFLOAT_IS_ZERO(x))
    {
        if (c == 0)
            return decfloat_zero(res, ctx);
        radix_integer_fit_limbs(&res->m, 1, DECIMAL_CTX_RADIX(ctx))[0] = c;
        res->m.size = 1;
        fmpz_set_si(&res->exp, k);
        return _decfloat_finalize(res, ctx);
    }

    n = FLINT_ABS(x->m.size);
    negative = x->m.size < 0;

    if (COEFF_IS_MPZ(x->exp))
        return GR_UNABLE;

    lo = FLINT_MIN(x->exp, k);
    hi = FLINT_MAX(x->exp + n, k + 1);
    span = hi - lo;
    if ((double) span * e > (double) DECIMAL_CONV_DIGITS_LIMIT)
        return GR_UNABLE;

    TMP_START;
    t = TMP_ALLOC(span * sizeof(ulong));
    flint_mpn_zero(t, span);
    flint_mpn_copyi(t + (x->exp - lo), x->m.d, n);
    t[k - lo] = c;

    /* strip zero limbs at both ends */
    for (a = 0; a < span && t[a] == 0; a++) ;
    for (b = span; b > a && t[b - 1] == 0; b--) ;

    if (a == b)
    {
        TMP_END;
        return decfloat_zero(res, ctx);
    }

    d = radix_integer_fit_limbs(&res->m, b - a, DECIMAL_CTX_RADIX(ctx));
    for (j = a; j < b; j++)
        d[j - a] = t[j];
    res->m.size = negative ? -(b - a) : (b - a);
    fmpz_set_si(&res->exp, lo + a);
    TMP_END;

    return _decfloat_finalize(res, ctx);
}

int
decfloat_set_digit_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);
    slong q = (k >= 0) ? k / e : -((-k + e - 1) / e);
    ulong p, old;

    if (c > 9)
        return GR_DOMAIN;

    p = n_pow(10, k - q * e);
    old = decfloat_get_limb_si(x, q, ctx);
    return decfloat_set_limb_si(res, x, q, old - ((old / p) % 10) * p + c * p, ctx);
}

int
decfloat_get_radix_integer_Bexp_fmpz(radix_integer_t m, fmpz_t exp, const decfloat_t x, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (!DECFLOAT_IS_ZERO(x))
            return DECFLOAT_IS_NAN(x) ? GR_UNABLE : GR_DOMAIN;
        radix_integer_zero(m, radix);
        fmpz_zero(exp);
        return GR_SUCCESS;
    }

    radix_integer_set(m, &x->m, radix);
    fmpz_set(exp, &x->exp);
    return GR_SUCCESS;
}

int
decfloat_get_radix_integer(radix_integer_t res, const decfloat_t x, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong n, s;
    nn_ptr d;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (!DECFLOAT_IS_ZERO(x))
            return DECFLOAT_IS_NAN(x) ? GR_UNABLE : GR_DOMAIN;
        radix_integer_zero(res, radix);
        return GR_SUCCESS;
    }

    if (fmpz_sgn(&x->exp) < 0)
        return GR_DOMAIN;

    n = FLINT_ABS(x->m.size);
    if (COEFF_IS_MPZ(x->exp) || ((double) x->exp + n) * DECIMAL_CTX_E(ctx) > (double) DECIMAL_CONV_DIGITS_LIMIT)
        return GR_UNABLE;

    s = x->exp;
    if (res == &x->m)
    {
        radix_integer_lshift_limbs(res, res, s, radix);
        return GR_SUCCESS;
    }
    d = radix_integer_fit_limbs(res, n + s, radix);
    flint_mpn_zero(d, s);
    flint_mpn_copyi(d + s, x->m.d, n);
    res->size = (x->m.size < 0) ? -(n + s) : (n + s);
    return GR_SUCCESS;
}

int
decfloat_set_radix_integer_Bexp_fmpz(decfloat_t res, const radix_integer_t m, const fmpz_t exp, gr_ctx_t ctx)
{
    slong n = FLINT_ABS(m->size);

    if (n == 0)
        return decfloat_zero(res, ctx);

    return _decfloat_set_round_limbs(res, m->d, n, m->size < 0, exp, 0,
        DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

int
decfloat_set_radix_integer(decfloat_t res, const radix_integer_t m, gr_ctx_t ctx)
{
    fmpz zero = 0;
    return decfloat_set_radix_integer_Bexp_fmpz(res, m, &zero, ctx);
}

void
decfloat_get_sci_exp(fmpz_t E, const decfloat_t x, gr_ctx_t ctx)
{
    FLINT_ASSERT(x->m.size != 0);
    fmpz_mul_ui(E, &x->exp, DECIMAL_CTX_E(ctx));
    fmpz_add_si(E, E, _decfloat_mant_digits(x, ctx) - 1);
}

int
decfloat_get_sci_exp_si(slong * E, const decfloat_t x, gr_ctx_t ctx)
{
    slong e = DECIMAL_CTX_E(ctx);

    FLINT_ASSERT(x->m.size != 0);

    if (COEFF_IS_MPZ(x->exp) || x->exp > WORD_MAX / (2 * e) || x->exp < -(WORD_MAX / (2 * e)))
        return 0;

    *E = e * x->exp + _decfloat_mant_digits(x, ctx) - 1;
    return 1;
}

/* ------------------------------------------------------------------------- */
/*    Random generation                                                      */
/* ------------------------------------------------------------------------- */

int
decfloat_randtest_special(decfloat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong prec = DECIMAL_CTX_PREC(ctx);
    slong maxlimbs;
    int status;

    if (prec == DECIMAL_PREC_EXACT)
        maxlimbs = 1 + n_randint(state, 4);
    else
        maxlimbs = 1 + n_randint(state, 2 + (prec + radix->exp - 1) / radix->exp);

    switch (n_randint(state, 16))
    {
        case 0:
            return decfloat_zero(res, ctx);
        case 1:
            status = decfloat_one(res, ctx);
            break;
        case 2:
            status = decfloat_neg_one(res, ctx);
            break;
        case 3:
            if (DECIMAL_CTX_ALLOW_INF(ctx))
                return n_randint(state, 2) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
            if (DECIMAL_CTX_ALLOW_NAN(ctx))
                return decfloat_nan(res, ctx);
            return decfloat_zero(res, ctx);
        case 4:
            status = decfloat_set_round_si(res, (slong) n_randtest(state), DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
            break;
        case 5:
            status = decfloat_set_round_si(res, (slong) n_randint(state, 21) - 10, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
            break;
        default:
            radix_integer_randtest_limbs(&res->m, state, maxlimbs, radix);

            if (res->m.size == 0)
                return decfloat_zero(res, ctx);

            /* strip low zero limbs */
            {
                slong n = FLINT_ABS(res->m.size), off = 0;
                while (res->m.d[off] == 0)
                    off++;
                if (off > 0)
                {
                    flint_mpn_copyi(res->m.d, res->m.d + off, n - off);
                    res->m.size = (res->m.size > 0) ? n - off : -(n - off);
                }
            }

            switch (n_randint(state, 8))
            {
                case 0:
                    /* huge exponents (avoiding the range where exact
                       conversions are feasible but very expensive) */
                    fmpz_randtest(&res->exp, state, 100);
                    if (fmpz_bits(&res->exp) > 8)
                        fmpz_mul_2exp(&res->exp, &res->exp, 40);
                    break;
                case 1:
                    fmpz_set_si(&res->exp, (slong) n_randint(state, 2 * maxlimbs + 2) - maxlimbs - 1);
                    break;
                default:
                    fmpz_set_si(&res->exp, (slong) n_randint(state, 2 * maxlimbs) - maxlimbs);
                    break;
            }

            status = GR_SUCCESS;
    }

    if (status != GR_SUCCESS)
        return decfloat_zero(res, ctx);

    return GR_SUCCESS;
}

int
decfloat_randtest(decfloat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;

    status = decfloat_randtest_special(res, state, ctx);

    if (status == GR_SUCCESS && !DECFLOAT_IS_SPECIAL(res))
    {
        int rnd = n_randint(state, DECIMAL_RND_NUM);
        status = decfloat_set_round(res, res, DECIMAL_CTX_PREC(ctx), rnd, ctx);
    }

    if (status != GR_SUCCESS)
        return decfloat_zero(res, ctx);

    return GR_SUCCESS;
}

POP_OPTIONS

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

static void
_set_rounding_exact(decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    if (info != NULL)
    {
        info->inexact = 0;
        info->increased = 0;
        info->underflow = 0;
        info->overflow = 0;
    }
    if (err != NULL)
        _decmag_zero(err, ctx);
}

static int
_decfloat_div_special(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    int sx, sy;

    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return decfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(y))
    {
        if (DECFLOAT_IS_ZERO(x))
            return DECIMAL_CTX_ALLOW_NAN(ctx) ? decfloat_nan(res, ctx) : GR_DOMAIN;
        if (!DECIMAL_CTX_ALLOW_INF(ctx))
            return GR_DOMAIN;
        return (_decfloat_sgn(x, ctx) > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
    }

    if (DECFLOAT_IS_INF(x) && DECFLOAT_IS_INF(y))
        return decfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(x) || DECFLOAT_IS_INF(y))
    {
        DECFLOAT_CHECK_OPERAND(x, ctx);
        DECFLOAT_CHECK_OPERAND(y, ctx);
        return decfloat_zero(res, ctx);
    }

    /* x is inf, y finite nonzero */
    sx = _decfloat_sgn(x, ctx);
    sy = _decfloat_sgn(y, ctx);
    return (sx * sy > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
}

/* Exact division in Z[1/10]: succeeds iff the reduced denominator has
   only the prime factors 2 and 5. */
static int
_decfloat_div_exact(decfloat_t res, const decfloat_t x, const decfloat_t y, int rnd, gr_ctx_t ctx)
{
    fmpz_t a, b, g, t, five, two;
    slong v2, v5, k;
    int status;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(g);
    fmpz_init(t);
    fmpz_init_set_ui(five, 5);
    fmpz_init_set_ui(two, 2);

    radix_integer_get_fmpz(a, &x->m, DECIMAL_CTX_RADIX(ctx));
    radix_integer_get_fmpz(b, &y->m, DECIMAL_CTX_RADIX(ctx));
    if (fmpz_sgn(b) < 0)
    {
        fmpz_neg(a, a);
        fmpz_neg(b, b);
    }

    fmpz_gcd(g, a, b);
    fmpz_divexact(a, a, g);
    fmpz_divexact(b, b, g);

    v2 = fmpz_remove(b, b, two);
    v5 = fmpz_remove(b, b, five);

    if (!fmpz_is_one(b))
    {
        status = GR_UNABLE;
    }
    else
    {
        /* a / (2^v2 5^v5) = a * 2^(k-v2) * 5^(k-v5) / 10^k, k = max(v2, v5) */
        k = FLINT_MAX(v2, v5);
        fmpz_ui_pow_ui(t, 2, k - v2);
        fmpz_mul(a, a, t);
        fmpz_ui_pow_ui(t, 5, k - v5);
        fmpz_mul(a, a, t);

        /* exponent: -k + e * (vx - vy) */
        fmpz_sub(t, &x->exp, &y->exp);
        fmpz_mul_ui(t, t, DECIMAL_CTX_E(ctx));
        fmpz_sub_ui(t, t, k);

        status = decfloat_set_round_fmpz_10exp_fmpz(res, a, t, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | (rnd & DECIMAL_RND_NOLIMITS), ctx);
    }

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(t);
    fmpz_clear(five);
    fmpz_clear(two);

    return status;
}

int
_decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong xn, yn, k, an, qn, i;
    int negative, status, sticky;
    fmpz_t exp;
    nn_ptr a, q, r;
    TMP_INIT;

    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y))
    {
        _set_rounding_exact(info, err, ctx);
        return _decfloat_div_special(res, x, y, ctx);
    }

    if (prec == DECIMAL_PREC_EXACT)
    {
        _set_rounding_exact(info, err, ctx);
        return _decfloat_div_exact(res, x, y, rnd, ctx);
    }

    xn = FLINT_ABS(x->m.size);
    yn = FLINT_ABS(y->m.size);
    negative = (x->m.size < 0) ^ (y->m.size < 0);

    /* Q = floor(Mx B^k / My) with at least prec + 1 digits */
    k = yn - xn + 1 + (prec + e - 1) / e;
    if (k < 0)
        k = 0;

    an = xn + k;
    qn = an - yn + 1;

    fmpz_init(exp);
    _fmpz_sub2_fast(exp, &x->exp, &y->exp, -k);

    TMP_START;
    a = TMP_ALLOC(sizeof(ulong) * (an + qn + 1 + yn));
    q = a + an;
    r = q + qn + 1;

    flint_mpn_zero(a, k);
    flint_mpn_copyi(a + k, x->m.d, xn);

    radix_divrem(q, r, a, an, y->m.d, yn, radix);

    sticky = 0;
    for (i = 0; i < yn; i++)
        sticky |= (r[i] != 0);

    status = _decfloat_set_round_limbs(res, q, qn, negative, exp, sticky, prec, rnd, info, err, ctx);

    TMP_END;
    fmpz_clear(exp);
    return status;
}

int
decfloat_div_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_div(res, x, y, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    return _decfloat_div(res, x, y, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

int
decfloat_inv_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t one;
    int status;

    decfloat_init(one, ctx);
    radix_integer_fit_limbs(&one->m, 1, DECIMAL_CTX_RADIX(ctx))[0] = 1;
    one->m.size = 1;
    status = _decfloat_div(res, one, x, prec, rnd, NULL, NULL, ctx);
    decfloat_clear(one, ctx);
    return status;
}

int
decfloat_inv(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_inv_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

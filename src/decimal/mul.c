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
#include "gr_generic.h"

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
_decfloat_mul_special(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    int sx, sy;

    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
        return decfloat_nan(res, ctx);

    if (DECFLOAT_IS_ZERO(x) || DECFLOAT_IS_ZERO(y))
    {
        if (DECFLOAT_IS_INF(x) || DECFLOAT_IS_INF(y))
            return decfloat_nan(res, ctx);
        return decfloat_zero(res, ctx);
    }

    sx = _decfloat_sgn(x, ctx);
    sy = _decfloat_sgn(y, ctx);

    return (sx * sy > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
}

/*
    Attempt to compute the rounded product from a truncated (high) product.
    The high product r approximates T = floor(xy / B^lo) with
    r <= T <= r + E, E = min(xn, yn, lo) * B. We round both r and r + E;
    if the two roundings agree, that is the correctly rounded product (by
    monotonicity of rounding). Returns 1 on success.
*/
static int
_decfloat_mul_high(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    ulong B = LIMB_RADIX(radix);
    slong xn = FLINT_ABS(x->m.size);
    slong yn = FLINT_ABS(y->m.size);
    slong pn = xn + yn;
    slong K, lo, hn, n1, n2, off1, off2, newn1, newn2, i;
    int negative = (x->m.size < 0) ^ (y->m.size < 0);
    nn_ptr t1, t2;
    ulong E, cy;
    fmpz_t exp;
    decimal_rounding_info info1, info2;
    decmag_t err1, err2;
    int ok, status;
    TMP_INIT;

    K = (prec + e - 1) / e + 2;
    lo = pn - K - 2;

    if (lo < 1)
        return 0;

    hn = pn - lo;

    /* the error bound E = min(xn, yn, lo) * B must fit in the two guard
       limbs; this can only fail for tiny limb radices */
    E = FLINT_MIN(FLINT_MIN(xn, yn), lo);
    if (E >= B)
        return 0;

    TMP_START;
    t1 = TMP_ALLOC(sizeof(ulong) * (hn + 1));
    t2 = TMP_ALLOC(sizeof(ulong) * (hn + 1));

    radix_mulmid(t1, x->m.d, xn, y->m.d, yn, lo, pn, radix);

    /* t2 = t1 + E where E = min(xn, yn, lo) * B */
    flint_mpn_copyi(t2, t1, hn);
    cy = radix_add(t2 + 1, t2 + 1, hn - 1, &E, 1, radix);
    t2[hn] = cy;

    n1 = hn;
    while (n1 > 0 && t1[n1 - 1] == 0)
        n1--;
    n2 = hn + 1;
    while (n2 > 0 && t2[n2 - 1] == 0)
        n2--;

    if (n1 == 0)
    {
        TMP_END;
        return 0;
    }

    fmpz_init(exp);
    _fmpz_add2_fast(exp, &x->exp, &y->exp, lo);

    _decmag_init(err1, ctx);
    _decmag_init(err2, ctx);

    off1 = _decimal_round_mantissa(t1, n1, negative, 0, prec, rnd, &newn1, &info1, err != NULL ? err1 : NULL, exp, ctx);
    off2 = _decimal_round_mantissa(t2, n2, negative, 0, prec, rnd, &newn2, &info2, err != NULL ? err2 : NULL, exp, ctx);

    ok = (off1 == off2 && newn1 == newn2);
    for (i = 0; ok && i < newn1; i++)
        ok = (t1[off1 + i] == t2[off2 + i]);

    if (ok)
    {
        nn_ptr rd = radix_integer_fit_limbs(&res->m, newn1, radix);
        flint_mpn_copyi(rd, t1 + off1, newn1);
        res->m.size = negative ? -newn1 : newn1;
        _fmpz_add_fast(&res->exp, exp, off1);

        if (info != NULL)
        {
            info->inexact = info1.inexact | info2.inexact;
            info->increased = info1.increased;
        }

        if (err != NULL)
            _decmag_max(err, err1, err2, ctx);

        if (rnd & DECIMAL_RND_NOLIMITS)
        {
            if (info != NULL)
            {
                info->underflow = 0;
                info->overflow = 0;
            }
            status = GR_SUCCESS;
        }
        else
            status = _decfloat_finalize_info(res, info, ctx);

        if (status != GR_SUCCESS)
            ok = -1;
    }

    _decmag_clear(err1, ctx);
    _decmag_clear(err2, ctx);
    fmpz_clear(exp);
    TMP_END;

    return ok;
}

int
_decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong xn, yn, pn;
    int negative, status;
    fmpz_t exp;
    nn_ptr t;

    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y))
    {
        _set_rounding_exact(info, err, ctx);
        return _decfloat_mul_special(res, x, y, ctx);
    }

    xn = FLINT_ABS(x->m.size);
    yn = FLINT_ABS(y->m.size);
    negative = (x->m.size < 0) ^ (y->m.size < 0);
    pn = xn + yn;

    /* single-limb product */
    if (pn == 2)
    {
        ulong q[2], t2[3];

        umul_ppmm(q[1], q[0], x->m.d[0], y->m.d[0]);
        t2[0] = flint_mpn_divrem_1_preinv(q, q, 2, radix->B.n, radix->B.ninv, radix->B.norm);
        t2[1] = q[0];   /* the product is < B^2, so the quotient fits in a limb */

        fmpz_init(exp);
        _fmpz_add2_fast(exp, &x->exp, &y->exp, 0);
        status = _decfloat_set_round_limbs(res, t2, 2, negative, exp, 0, prec, rnd, info, err, ctx);
        fmpz_clear(exp);
        return status;
    }

    /* truncated product when only the top limbs are needed */
    if (prec != DECIMAL_PREC_EXACT)
    {
        int r = _decfloat_mul_high(res, x, y, prec, rnd, info, err, ctx);
        if (r == 1)
            return GR_SUCCESS;
        if (r == -1)
            return GR_UNABLE;
    }

    /* refuse to materialize huge exact products */
    if (prec == DECIMAL_PREC_EXACT && pn > DECIMAL_CONV_DIGITS_LIMIT / radix->exp)
        return GR_UNABLE;

    fmpz_init(exp);
    _fmpz_add2_fast(exp, &x->exp, &y->exp, 0);

    if (res != x && res != y)
    {
        t = radix_integer_fit_limbs(&res->m, pn + 1, radix);

        if (x == y)
            radix_sqr(t, x->m.d, xn, radix);
        else
            radix_mul(t, x->m.d, xn, y->m.d, yn, radix);

        status = _decfloat_set_round_limbs(res, t, pn, negative, exp, 0, prec, rnd, info, err, ctx);
    }
    else
    {
        TMP_INIT;
        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * (pn + 1));

        if (x == y)
            radix_sqr(t, x->m.d, xn, radix);
        else
            radix_mul(t, x->m.d, xn, y->m.d, yn, radix);

        status = _decfloat_set_round_limbs(res, t, pn, negative, exp, 0, prec, rnd, info, err, ctx);
        TMP_END;
    }

    fmpz_clear(exp);
    return status;
}

int
decfloat_mul_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_mul(res, x, y, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_sqr_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_mul(res, x, x, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    return _decfloat_mul(res, x, y, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

int
decfloat_sqr(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return _decfloat_mul(res, x, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

DECIMAL_DRIVER int
_decfloat_mul_scalar(decfloat_t res, const decfloat_t x, const void * y, int type, int div, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = _decfloat_set_scalar_exact(t, y, type, ctx);
    if (status == GR_SUCCESS)
        status = (div ? _decfloat_div : _decfloat_mul)(res, x, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int decfloat_mul_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, &y, DECIMAL_SCALAR_UI, 0, ctx); }
int decfloat_mul_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, &y, DECIMAL_SCALAR_SI, 0, ctx); }
int decfloat_mul_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, y, DECIMAL_SCALAR_FMPZ, 0, ctx); }
int decfloat_div_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, &y, DECIMAL_SCALAR_UI, 1, ctx); }
int decfloat_div_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, &y, DECIMAL_SCALAR_SI, 1, ctx); }
int decfloat_div_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx) { return _decfloat_mul_scalar(res, x, y, DECIMAL_SCALAR_FMPZ, 1, ctx); }

/* binary powering with rounded multiplications (exact contexts and
   special values); does not depend on the method table of ctx, so that
   it can be called with a complex context */
static int
_decfloat_pow_fmpz_binexp(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
{
    decfloat_t t;
    slong i, bits;
    int status = GR_SUCCESS;

    if (fmpz_is_zero(y))
    {
        DECFLOAT_CHECK_OPERAND(x, ctx);
        return decfloat_one(res, ctx);
    }

    if (fmpz_is_one(y))
        return decfloat_set(res, x, ctx);

    /* an exact power of a number with d > 1 significant digits (or a
       single digit other than 1) has about d |y| digits */
    if (DECIMAL_CTX_IS_EXACT(ctx) && !DECFLOAT_IS_SPECIAL(x))
    {
        slong d = decfloat_digits(x, ctx);
        ulong v = x->m.d[0];
        while (v % 10 == 0)
            v /= 10;
        if (d > 1 || v != 1)
            if (fmpz_bits(y) > 40 || (double) FLINT_MAX(d - 1, 1) * fmpz_get_d(y) * (fmpz_sgn(y) < 0 ? -1 : 1) > (double) DECIMAL_CONV_DIGITS_LIMIT)
                return GR_UNABLE;
    }

    decfloat_init(t, ctx);
    status = decfloat_set(t, x, ctx);

    /* bits of |y| */
    {
        fmpz_t a;
        fmpz_init(a);
        fmpz_abs(a, y);
        bits = fmpz_bits(a);
        for (i = bits - 2; i >= 0 && status == GR_SUCCESS; i--)
        {
            status = decfloat_mul(t, t, t, ctx);
            if (status == GR_SUCCESS && fmpz_tstbit(a, i))
                status = decfloat_mul(t, t, x, ctx);
        }
        fmpz_clear(a);
    }

    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(y) < 0)
            status = decfloat_inv(res, t, ctx);
        else
            decfloat_swap(res, t, ctx);
    }

    decfloat_clear(t, ctx);
    return status;
}

int
decfloat_pow_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
{
    int r;

    if (DECFLOAT_IS_SPECIAL(x) || fmpz_is_zero(y))
        return _decfloat_pow_fmpz_binexp(res, x, y, ctx);

    /* exact power followed by a single rounding when feasible */
    r = _decfloat_pow_int_exact(res, x, y, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
    if (r == 1)
        return GR_SUCCESS;
    if (r == -1)
        return GR_UNABLE;

    if (DECIMAL_CTX_IS_EXACT(ctx))
        return _decfloat_pow_fmpz_binexp(res, x, y, ctx);

    /* correctly rounded via arb */
    {
        decfloat_t t;
        int status;

        if (_decfloat_sgn(x, ctx) < 0)
        {
            /* (-x)^n = (-1)^n x^n */
            decfloat_t u;
            decfloat_init(u, ctx);
            GR_MUST_SUCCEED(decfloat_neg_round(u, x, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
            if (fmpz_is_odd(y))
            {
                /* -(x^n) rounded correctly: use the mirrored rounding mode */
                int saved = DECIMAL_CTX_RND(ctx);
                DECIMAL_CTX_RND(ctx) = DECIMAL_RND_NEGATE(saved);
                status = decfloat_pow_fmpz(res, u, y, ctx);
                DECIMAL_CTX_RND(ctx) = saved;
                if (status == GR_SUCCESS)
                    status = decfloat_neg_round(res, res, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx);
            }
            else
                status = decfloat_pow_fmpz(res, u, y, ctx);
            decfloat_clear(u, ctx);
            return status;
        }

        decfloat_init(t, ctx);
        GR_MUST_SUCCEED(decfloat_set_round_fmpz(t, y, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN | DECIMAL_RND_NOLIMITS, ctx));
        status = _decfloat_pow_ziv(res, x, t, ctx);
        decfloat_clear(t, ctx);
        return status;
    }
}

int
decfloat_pow_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_ui(t, y);
    status = decfloat_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

int
decfloat_pow_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init_set_si(t, y);
    status = decfloat_pow_fmpz(res, x, t, ctx);
    fmpz_clear(t);
    return status;
}

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

int
_decfloat_mul_10exp(decfloat_t res, const decfloat_t x, const fmpz_t t, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    fmpz_t q, exp;
    ulong r;
    slong n;
    int negative, status;

    if (DECFLOAT_IS_SPECIAL(x) || fmpz_is_zero(t))
        return decfloat_set_round_info(res, x, prec, rnd, info, err, ctx);

    fmpz_init(q);
    fmpz_init(exp);

    if (!COEFF_IS_MPZ(*t) && *t >= 0 && *t < e)
    {
        r = *t;
    }
    else
    {
        fmpz_t rr;
        fmpz_init(rr);
        fmpz_fdiv_q_ui(q, t, e);
        fmpz_mul_ui(rr, q, e);
        fmpz_sub(rr, t, rr);
        r = fmpz_get_ui(rr);
        fmpz_clear(rr);
    }

    fmpz_add(exp, &x->exp, q);

    n = FLINT_ABS(x->m.size);
    negative = x->m.size < 0;

    if (r == 0)
    {
        if (res != x)
        {
            nn_ptr d = radix_integer_fit_limbs(&res->m, n, radix);
            flint_mpn_copyi(d, x->m.d, n);
            res->m.size = x->m.size;
        }
        fmpz_swap(&res->exp, exp);
        status = decfloat_set_round_info(res, res, prec, rnd, info, err, ctx);
    }
    else
    {
        nn_ptr d;
        ulong cy;

        if (res == x)
            d = radix_integer_fit_limbs(&res->m, n + 2, radix);
        else
        {
            d = radix_integer_fit_limbs(&res->m, n + 2, radix);
        }

        cy = radix_lshift_digits(d, x->m.d, n, r, radix);
        d[n] = cy;

        status = _decfloat_set_round_limbs(res, d, n + 1, negative, exp, 0, prec, rnd, info, err, ctx);
    }

    fmpz_clear(q);
    fmpz_clear(exp);
    return status;
}

int
decfloat_mul_10exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t t, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_mul_10exp(res, x, t, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_mul_10exp_si_round(decfloat_t res, const decfloat_t x, slong t, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz_t tt;
    int status;
    fmpz_init_set_si(tt, t);
    status = _decfloat_mul_10exp(res, x, tt, prec, rnd, NULL, NULL, ctx);
    fmpz_clear(tt);
    return status;
}

int
decfloat_mul_10exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t t, gr_ctx_t ctx)
{
    return _decfloat_mul_10exp(res, x, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

int
decfloat_mul_10exp_si(decfloat_t res, const decfloat_t x, slong t, gr_ctx_t ctx)
{
    return decfloat_mul_10exp_si_round(res, x, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_mul_2exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t t, slong prec, int rnd, gr_ctx_t ctx)
{
    fmpz_t m, e, p;
    int status;

    if (DECFLOAT_IS_SPECIAL(x) || fmpz_is_zero(t))
        return decfloat_set_round(res, x, prec, rnd, ctx);

    if (prec == DECIMAL_PREC_EXACT)
    {
        /* the exact result has about digits(x) + 0.7 |t| digits */
        if (fmpz_bits(t) > 40 || (double) decfloat_digits(x, ctx) + 0.7 * fmpz_get_d(t) * (fmpz_sgn(t) < 0 ? -1 : 1) > (double) DECIMAL_CONV_DIGITS_LIMIT)
            return GR_UNABLE;
    }

    fmpz_init(m);
    fmpz_init(e);
    fmpz_init(p);

    /* x = m * 10^e exactly */
    GR_MUST_SUCCEED(decfloat_get_fmpz_10exp_fmpz(m, e, x, ctx));

    if (prec != DECIMAL_PREC_EXACT && fmpz_bits(t) > 16)
    {
        /* large |t|: correctly round m 2^t with the binary conversion
           (without the exponent limits), then shift by 10^e exactly and
           apply the limits (overflow and underflow are decided by the
           unconstrained rounded value) */
        status = decfloat_set_round_fmpz_2exp_fmpz(res, m, t, prec, rnd | DECIMAL_RND_NOLIMITS, ctx);
        if (status == GR_SUCCESS)
            status = decfloat_mul_10exp_fmpz_round(res, res, e, prec, rnd, ctx);
        fmpz_clear(m);
        fmpz_clear(e);
        fmpz_clear(p);
        return status;
    }

    if (fmpz_sgn(t) > 0)
    {
        fmpz_mul_2exp(m, m, fmpz_get_ui(t));
    }
    else
    {
        /* 2^-k = 5^k 10^-k */
        ulong k = -fmpz_get_si(t);
        fmpz_ui_pow_ui(p, 5, k);
        fmpz_mul(m, m, p);
        fmpz_sub_ui(e, e, k);
    }

    status = decfloat_set_round_fmpz_10exp_fmpz(res, m, e, prec, rnd, ctx);

    fmpz_clear(m);
    fmpz_clear(e);
    fmpz_clear(p);
    return status;
}

int
decfloat_mul_2exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t t, gr_ctx_t ctx)
{
    return decfloat_mul_2exp_fmpz_round(res, x, t, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_mul_2exp_si(decfloat_t res, const decfloat_t x, slong t, gr_ctx_t ctx)
{
    fmpz_t tt;
    int status;
    fmpz_init_set_si(tt, t);
    status = decfloat_mul_2exp_fmpz(res, x, tt, ctx);
    fmpz_clear(tt);
    return status;
}

/* ------------------------------------------------------------------------- */
/*    Rounding to integers                                                   */
/* ------------------------------------------------------------------------- */

int
_decfloat_round_to_int(decfloat_t res, const decfloat_t x, int int_rnd, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong n, k, D_int;
    int negative, status;

    if (DECFLOAT_IS_SPECIAL(x) || fmpz_sgn(&x->exp) >= 0)
        return decfloat_set_round_info(res, x, prec, rnd, info, err, ctx);

    n = FLINT_ABS(x->m.size);
    negative = x->m.size < 0;

    /* all limbs fractional: |x| < 1 */
    if (fmpz_cmp_si(&x->exp, -n) <= 0)
    {
        int result = 0;   /* -1, 0, 1 */

        switch (int_rnd)
        {
            case DECIMAL_RND_FLOOR:
                result = negative ? -1 : 0;
                break;
            case DECIMAL_RND_CEIL:
                result = negative ? 0 : 1;
                break;
            case DECIMAL_RND_DOWN:
                result = 0;
                break;
            case DECIMAL_RND_UP:
                result = negative ? -1 : 1;
                break;
            default:
                {
                    /* compare |x| with 1/2 */
                    int c;
                    if (fmpz_cmp_si(&x->exp, -n) < 0)
                        c = -1;
                    else
                    {
                        ulong half = LIMB_RADIX(radix) / 2;
                        slong i;
                        c = (x->m.d[n - 1] < half) ? -1 : (x->m.d[n - 1] > half ? 1 : 0);
                        if (c == 0)
                            for (i = 0; i < n - 1; i++)
                                if (x->m.d[i] != 0) { c = 1; break; }
                    }

                    if (c > 0 || (c == 0 && int_rnd == DECIMAL_RND_NEAR_AWAY))
                        result = negative ? -1 : 1;
                    else
                        result = 0;
                }
        }

        if (info != NULL)
        {
            info->inexact = 1;
            info->increased = (result != 0);
            info->underflow = 0;
            info->overflow = 0;
        }
        if (err != NULL)
            _decmag_one(err, ctx);

        if (result == 0)
            return decfloat_zero(res, ctx);

        return decfloat_set_round_si(res, result, prec, rnd, ctx);
    }

    k = -fmpz_get_si(&x->exp);    /* 1 <= k < n */
    D_int = (n - k - 1) * e + _radix_size_digits_1(x->m.d[n - 1], radix);

    {
        nn_ptr d = radix_integer_fit_limbs(&res->m, n + 1, radix);
        fmpz_t zero_exp;

        if (d != x->m.d)
            flint_mpn_copyi(d, x->m.d, n);

        fmpz_init(zero_exp);
        fmpz_set(zero_exp, &x->exp);

        status = _decfloat_set_round_limbs(res, d, n, negative, zero_exp, 0, D_int, int_rnd | DECIMAL_RND_NOLIMITS, info, err, ctx);
        fmpz_clear(zero_exp);
    }

    if (status == GR_SUCCESS)
        status = decfloat_set_round_info(res, res, prec, rnd, NULL, NULL, ctx);

    return status;
}

int
decfloat_floor_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_round_to_int(res, x, DECIMAL_RND_FLOOR, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_ceil_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_round_to_int(res, x, DECIMAL_RND_CEIL, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_trunc_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_round_to_int(res, x, DECIMAL_RND_DOWN, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_nint_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_round_to_int(res, x, DECIMAL_RND_NEAR, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_floor(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_floor_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_ceil(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_ceil_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_trunc(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_trunc_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

int
decfloat_nint(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_nint_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

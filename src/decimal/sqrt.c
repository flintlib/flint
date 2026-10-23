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

int
_decfloat_sqrt(decfloat_t res, const decfloat_t x, slong prec, int rnd,
    decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    slong xn, an, sn, k, i, shift;
    int status, sticky, odd;
    fmpz_t exp;
    nn_ptr a, s, r;
    TMP_INIT;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        _set_rounding_exact(info, err, ctx);
        if (DECFLOAT_IS_ZERO(x))
            return decfloat_zero(res, ctx);
        if (DECFLOAT_IS_POS_INF(x))
            return decfloat_pos_inf(res, ctx);
        if (DECFLOAT_IS_NEG_INF(x))
            return GR_DOMAIN;
        return decfloat_nan(res, ctx);
    }

    if (x->m.size < 0)
        return GR_DOMAIN;

    xn = x->m.size;
    odd = fmpz_is_odd(&x->exp);

    fmpz_init(exp);

    if (prec == DECIMAL_PREC_EXACT)
    {
        radix_integer_t t, u;
        int square;

        _set_rounding_exact(info, err, ctx);

        radix_integer_init(t, radix);
        radix_integer_init(u, radix);

        if (odd)
            radix_integer_lshift_limbs(t, &x->m, 1, radix);
        else
            radix_integer_set(t, &x->m, radix);

        square = radix_integer_sqrt(u, t, radix);

        if (square)
        {
            fmpz_sub_ui(exp, &x->exp, odd);
            fmpz_fdiv_q_2exp(exp, exp, 1);
            FLINT_SWAP(radix_integer_struct, res->m, *u);
            fmpz_swap(&res->exp, exp);
            /* the root has no low zero limbs if the input had none... but
               with odd shift it may; normalize */
            {
                slong n = res->m.size, off = 0;
                while (res->m.d[off] == 0)
                    off++;
                if (off > 0)
                {
                    flint_mpn_copyi(res->m.d, res->m.d + off, n - off);
                    res->m.size = n - off;
                    fmpz_add_ui(&res->exp, &res->exp, off);
                }
            }
            status = _decfloat_finalize(res, ctx);
        }
        else
        {
            status = GR_UNABLE;
        }

        radix_integer_clear(t, radix);
        radix_integer_clear(u, radix);
        fmpz_clear(exp);
        return status;
    }

    /* A = Mx * B^(odd + 2k), an = xn + odd + 2k >= 2 ceil(prec/e) + 3 */
    shift = 2 * ((prec + e - 1) / e) + 3 - (xn + odd);
    k = (shift <= 0) ? 0 : (shift + 1) / 2;
    an = xn + odd + 2 * k;
    sn = (an + 1) / 2;

    /* exponent: (v - odd)/2 - k */
    fmpz_sub_ui(exp, &x->exp, odd);
    fmpz_fdiv_q_2exp(exp, exp, 1);
    fmpz_sub_ui(exp, exp, k);

    TMP_START;
    a = TMP_ALLOC(sizeof(ulong) * (an + (sn + 1) + (sn + 2)));
    s = a + an;
    r = s + sn + 1;

    flint_mpn_zero(a, odd + 2 * k);
    flint_mpn_copyi(a + odd + 2 * k, x->m.d, xn);

    radix_sqrtrem(s, r, a, an, radix);

    sticky = 0;
    for (i = 0; i < sn + 1; i++)
        sticky |= (r[i] != 0);

    status = _decfloat_set_round_limbs(res, s, sn, 0, exp, sticky, prec, rnd, info, err, ctx);

    TMP_END;
    fmpz_clear(exp);
    return status;
}

int
decfloat_sqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_sqrt(res, x, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_sqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return _decfloat_sqrt(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

/* rsqrt(x) = 1 / sqrt(x): computed as a correctly rounded quotient
   after an exact-enough square root (Ziv style). */
int
decfloat_rsqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
{
    decfloat_t s, one, q;
    decimal_rounding_info info;
    slong wp;
    int status;

    if (DECFLOAT_IS_SPECIAL(x))
    {
        if (DECFLOAT_IS_ZERO(x))
            return DECIMAL_CTX_ALLOW_INF(ctx) ? decfloat_pos_inf(res, ctx) : GR_DOMAIN;
        if (DECFLOAT_IS_POS_INF(x))
        {
            DECFLOAT_CHECK_OPERAND(x, ctx);
            return decfloat_zero(res, ctx);
        }
        if (DECFLOAT_IS_NEG_INF(x))
            return GR_DOMAIN;
        return decfloat_nan(res, ctx);
    }

    if (x->m.size < 0)
        return GR_DOMAIN;

    if (prec == DECIMAL_PREC_EXACT)
    {
        decfloat_init(s, ctx);
        status = _decfloat_sqrt(s, x, DECIMAL_PREC_EXACT, rnd, NULL, NULL, ctx);
        if (status == GR_SUCCESS)
            status = decfloat_inv_round(res, s, DECIMAL_PREC_EXACT, rnd, ctx);
        decfloat_clear(s, ctx);
        return status;
    }

    decfloat_init(s, ctx);
    decfloat_init(one, ctx);
    decfloat_init(q, ctx);
    GR_MUST_SUCCEED(decfloat_one(one, ctx));

    /* 1/sqrt(x): compute sqrt to prec + guard digits, divide with the
       sticky flag from the square root, and verify that rounding is
       decided; otherwise increase the working precision. */
    status = GR_UNABLE;
    for (wp = prec + 10; wp < 100 * prec + 1000; wp *= 2)
    {
        status = _decfloat_sqrt(s, x, wp, DECIMAL_RND_DOWN, &info, NULL, ctx);
        if (status != GR_SUCCESS)
            break;

        if (!info.inexact)
        {
            /* exact square root: quotient is correctly rounded */
            status = _decfloat_div(res, one, s, prec, rnd, NULL, NULL, ctx);
            break;
        }
        else
        {
            /* s < sqrt(x) < s + ulp(s) at wp digits; 1/(s+ulp) < 1/sqrt(x) < 1/s.
               Round both endpoints at prec digits; if they agree, done. */
            decfloat_t s2, q2;
            decfloat_init(s2, ctx);
            decfloat_init(q2, ctx);

            status = _decfloat_div(q, one, s, prec, rnd, NULL, NULL, ctx);
            if (status == GR_SUCCESS)
            {
                decmag_t ulp;
                decfloat_t u;
                _decmag_init(ulp, ctx);
                decfloat_init(u, ctx);
                _decmag_set_ulp(ulp, s, wp, ctx);
                GR_MUST_SUCCEED(_decmag_get_decfloat(u, ulp, ctx));
                status = _decfloat_add(s2, s, u, 0, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, NULL, NULL, ctx);
                if (status == GR_SUCCESS)
                    status = _decfloat_div(q2, one, s2, prec, rnd, NULL, NULL, ctx);
                _decmag_clear(ulp, ctx);
                decfloat_clear(u, ctx);
            }

            if (status == GR_SUCCESS)
            {
                if (decfloat_equal(q, q2, ctx) == T_TRUE)
                {
                    decfloat_swap(res, q, ctx);
                    decfloat_clear(s2, ctx);
                    decfloat_clear(q2, ctx);
                    break;
                }
            }

            decfloat_clear(s2, ctx);
            decfloat_clear(q2, ctx);

            if (status != GR_SUCCESS)
                break;

            status = GR_UNABLE;
        }
    }

    decfloat_clear(s, ctx);
    decfloat_clear(one, ctx);
    decfloat_clear(q, ctx);
    return status;
}

int
decfloat_rsqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return decfloat_rsqrt_round(res, x, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), ctx);
}

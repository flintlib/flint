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
_decfloat_add_special(decfloat_t res, const decfloat_t x, const decfloat_t y, int negate_y,
    slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    if (DECFLOAT_IS_NAN(x) || DECFLOAT_IS_NAN(y))
    {
        _set_rounding_exact(info, err, ctx);
        return decfloat_nan(res, ctx);
    }

    if (DECFLOAT_IS_ZERO(x))
    {
        if (negate_y)
        {
            int status = decfloat_set_round_info(res, y, prec, DECIMAL_RND_NEGATE(rnd), info, err, ctx);
            if (status == GR_SUCCESS)
            {
                if (DECFLOAT_IS_POS_INF(res)) _decfloat_neg_inf(res);
                else if (DECFLOAT_IS_NEG_INF(res)) _decfloat_pos_inf(res);
                else res->m.size = -res->m.size;
            }
            return status;
        }
        return decfloat_set_round_info(res, y, prec, rnd, info, err, ctx);
    }

    if (DECFLOAT_IS_ZERO(y))
        return decfloat_set_round_info(res, x, prec, rnd, info, err, ctx);

    _set_rounding_exact(info, err, ctx);

    /* at least one infinity */
    {
        int sx = _decfloat_sgn(x, ctx);
        int sy = _decfloat_sgn(y, ctx);
        int xinf = DECFLOAT_IS_INF(x);
        int yinf = DECFLOAT_IS_INF(y);

        if (negate_y)
            sy = -sy;

        if (xinf && yinf)
        {
            if (sx != sy)
                return decfloat_nan(res, ctx);
            return (sx > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
        }

        if (xinf)
            return (sx > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
        else
            return (sy > 0) ? decfloat_pos_inf(res, ctx) : decfloat_neg_inf(res, ctx);
    }
}

int
_decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, int negate_y,
    slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
{
    const radix_struct * radix = DECIMAL_CTX_RADIX(ctx);
    slong e = radix->exp;
    ulong B = LIMB_RADIX(radix);
    const decfloat_struct * a;
    const decfloat_struct * b;
    slong an, bn, shift, K, span, offa, offb, tn;
    int sa, sb, negative, status, eps;
    nn_ptr t;
#define DECIMAL_ADD_STACK_LIMBS 8
    ulong tstack[DECIMAL_ADD_STACK_LIMBS];
    TMP_INIT;

    if (DECFLOAT_IS_SPECIAL(x) || DECFLOAT_IS_SPECIAL(y))
        return _decfloat_add_special(res, x, y, negate_y, prec, rnd, info, err, ctx);

    an = FLINT_ABS(x->m.size);
    bn = FLINT_ABS(y->m.size);
    sa = (x->m.size < 0);
    sb = (y->m.size < 0) ^ (negate_y != 0);

    /* top_x - top_y */
    shift = _fmpz_sub_small(&x->exp, &y->exp);
    if (shift < WORD_MAX / 4 && shift > -WORD_MAX / 4)
        shift += an - bn;

    if (shift >= 0)
    {
        a = x;
        b = y;
    }
    else
    {
        a = y;
        b = x;
        FLINT_SWAP(slong, an, bn);
        FLINT_SWAP(int, sa, sb);
        shift = -shift;
    }

    /* now a has the higher (or equal) top; shift = top_a - top_b >= 0 */

    if (prec == DECIMAL_PREC_EXACT)
    {
        K = DECIMAL_CONV_DIGITS_LIMIT / e;
        if (shift >= K)
            return GR_UNABLE;
    }
    else
    {
        K = (prec + e - 1) / e + 4;
    }

    /* b lies entirely (with a gap) below a's lowest limb and below the
       rounding horizon: a +/- tiny */
    if (shift > K && shift >= an + 1)
    {
        fmpz_t exp;
        slong pad;

        /* pad a with zero limbs so that the rounding position lies strictly
           above the sticky tail */
        pad = (prec + e - 1) / e + 2 - an;
        if (pad < 0)
            pad = 0;

        fmpz_init(exp);
        _fmpz_add_fast(exp, &a->exp, -pad);

        TMP_START;
        t = TMP_ALLOC(sizeof(ulong) * (an + pad + 1));
        flint_mpn_zero(t, pad);
        flint_mpn_copyi(t + pad, a->m.d, an);
        tn = an + pad;

        if (sa != sb)
        {
            /* a - tiny = (a - 1 unit) + (1 unit - tiny) */
            slong i = 0;
            while (t[i] == 0)
            {
                t[i] = B - 1;
                i++;
            }
            t[i] -= 1;
        }

        status = _decfloat_set_round_limbs(res, t, tn, sa, exp, 1, prec, rnd, info, err, ctx);
        TMP_END;
        fmpz_clear(exp);
        return status;
    }

    /* materialize */
    {
        /* lo = min(v_a, v_b) where v_a = top_a - an, v_b = top_b - bn = top_a - shift - bn */
        slong va_rel = 0;                 /* v_a relative to top_a */
        slong vb_rel = -shift - bn;
        slong lo_rel;
        const fmpz * lo_exp;

        va_rel = -an;

        if (va_rel <= vb_rel)
        {
            lo_rel = va_rel;
            lo_exp = &a->exp;
        }
        else
        {
            lo_rel = vb_rel;
            lo_exp = &b->exp;
        }

        span = -lo_rel;           /* limbs from lo to top_a */
        offa = va_rel - lo_rel;
        offb = vb_rel - lo_rel;

        TMP_START;
        if (span + 1 <= DECIMAL_ADD_STACK_LIMBS)
            t = tstack;
        else
            t = TMP_ALLOC(sizeof(ulong) * (span + 1));

        if (sa == sb)
        {
            /* L is the operand with offset 0, H the other one at offset oh */
            const decfloat_struct * L = (offa == 0) ? a : b;
            const decfloat_struct * H = (offa == 0) ? b : a;
            slong nL = (offa == 0) ? an : bn;
            slong nH = (offa == 0) ? bn : an;
            slong oh = (offa == 0) ? offb : offa;
            ulong cy;

            if (oh >= nL)
            {
                /* no overlap: concatenate */
                flint_mpn_copyi(t, L->m.d, nL);
                flint_mpn_zero(t + nL, oh - nL);
                flint_mpn_copyi(t + oh, H->m.d, nH);
                cy = 0;
            }
            else
            {
                flint_mpn_copyi(t, L->m.d, oh);
                if (nL - oh >= nH)
                    cy = radix_add(t + oh, L->m.d + oh, nL - oh, H->m.d, nH, radix);
                else
                    cy = radix_add(t + oh, H->m.d, nH, L->m.d + oh, nL - oh, radix);
            }

            t[span] = cy;
            tn = span + 1;
            negative = sa;
            eps = 0;
        }
        else
        {
            int cmp = _decfloat_cmpabs_finite(a, b);
            const decfloat_struct * big;
            const decfloat_struct * small;
            slong bign, smalln, bigoff, smalloff;

            if (cmp == 0)
            {
                TMP_END;
                _set_rounding_exact(info, err, ctx);
                return decfloat_zero(res, ctx);
            }

            if (cmp > 0)
            {
                big = a; small = b; bign = an; smalln = bn; bigoff = offa; smalloff = offb;
                negative = sa;
            }
            else
            {
                big = b; small = a; bign = bn; smalln = an; bigoff = offb; smalloff = offa;
                negative = sb;
            }

            flint_mpn_zero(t, bigoff);
            flint_mpn_copyi(t + bigoff, big->m.d, bign);
            flint_mpn_zero(t + bigoff + bign, span - bigoff - bign);
            radix_sub(t + smalloff, t + smalloff, span - smalloff, small->m.d, smalln, radix);
            tn = span;
            eps = 0;
        }

        status = _decfloat_set_round_limbs(res, t, tn, negative, lo_exp, eps, prec, rnd, info, err, ctx);
        TMP_END;
        return status;
    }
}

int
decfloat_add_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_add(res, x, y, 0, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_sub_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
{
    return _decfloat_add(res, x, y, 1, prec, rnd, NULL, NULL, ctx);
}

int
decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    return _decfloat_add(res, x, y, 0, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

int
decfloat_sub(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
{
    return _decfloat_add(res, x, y, 1, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

DECIMAL_DRIVER int
_decfloat_add_scalar(decfloat_t res, const decfloat_t x, const void * y, int type, int negate, gr_ctx_t ctx)
{
    decfloat_t t;
    int status;
    decfloat_init(t, ctx);
    status = _decfloat_set_scalar_exact(t, y, type, ctx);
    if (status == GR_SUCCESS)
        status = _decfloat_add(res, x, t, negate, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
    decfloat_clear(t, ctx);
    return status;
}

int decfloat_add_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, &y, DECIMAL_SCALAR_UI, 0, ctx); }
int decfloat_add_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, &y, DECIMAL_SCALAR_SI, 0, ctx); }
int decfloat_add_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, y, DECIMAL_SCALAR_FMPZ, 0, ctx); }
int decfloat_sub_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, &y, DECIMAL_SCALAR_UI, 1, ctx); }
int decfloat_sub_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, &y, DECIMAL_SCALAR_SI, 1, ctx); }
int decfloat_sub_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx) { return _decfloat_add_scalar(res, x, y, DECIMAL_SCALAR_FMPZ, 1, ctx); }

int
decfloat_mul_two(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
{
    return _decfloat_add(res, x, x, 0, DECIMAL_CTX_PREC(ctx), DECIMAL_CTX_RND(ctx), NULL, NULL, ctx);
}

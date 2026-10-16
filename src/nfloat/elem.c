/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "elem.h"

/* Result assembly ***********************************************************/

/* res = (-1)^sgnbit (y, len) 2^e truncated, (y, len) nonzero */
static int
_nfloat_set_mpn_trunc(nfloat_ptr res, nn_srcptr y, slong len, slong e, int sgnbit, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    slong i, exp;
    unsigned int c;
    nn_ptr r = NFLOAT_D(res);

    while (y[len - 1] == 0)
        len--;

    c = flint_clz(y[len - 1]);
    exp = e + FLINT_BITS * len - c;

    if (len > n)
    {
        y += len - n;

        if (c == 0)
        {
            for (i = 0; i < n; i++)
                r[i] = y[i];
        }
        else
        {
            for (i = 0; i < n; i++)
                r[i] = (y[i] << c) | (y[i - 1] >> (FLINT_BITS - c));
        }
    }
    else
    {
        slong k = n - len;

        for (i = 0; i < k; i++)
            r[i] = 0;

        if (c == 0)
        {
            for (i = 0; i < len; i++)
                r[k + i] = y[i];
        }
        else
        {
            r[k] = y[0] << c;
            for (i = 1; i < len; i++)
                r[k + i] = (y[i] << c) | (y[i - 1] >> (FLINT_BITS - c));
        }
    }

    NFLOAT_EXP(res) = exp;
    NFLOAT_SGNBIT(res) = sgnbit;
    NFLOAT_HANDLE_UNDERFLOW_OVERFLOW(res, ctx);
    return GR_SUCCESS;
}

int
_nfloat_set_mpn_err(nfloat_ptr res, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, gr_ctx_t ctx)
{
    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        return _nfloat_set_mpn_trunc(res, y, len, e, sgnbit, ctx);
    else
    {
        if (nfloat_should_round_up(sgnbit, ctx))
        {
            /* the endpoint away from zero, rounded away from zero */
            y[len] = mpn_add_1(y, y, len, ehi);
            len++;
        }
        else if (mpn_sub_1(y, y, len, elo))
        {
            /* the endpoint towards zero has the opposite sign: its
               magnitude elo - Y, rounded away from zero (which is the
               rounding direction for the flipped sign) */
            mpn_neg(y, y, len);
            sgnbit ^= 1;
        }
    }

    while (len > 0 && y[len - 1] == 0)
        len--;

    if (len == 0)
        return nfloat_zero(res, ctx);

    return _nfloat_set_mpn_2exp(res, y, len, e + FLINT_BITS * len, sgnbit, ctx);
}

int
_nfloat_one_plus_tiny(nfloat_ptr res, int sgnbit, int dsgn, slong e, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);

    FLINT_ASSERT(e <= -FLINT_BITS * n);
    (void) e;

    NFLOAT_SGNBIT(res) = sgnbit;

    if (NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
    {
        int up = nfloat_should_round_up(sgnbit, ctx);

        if (up && !dsgn)
        {
            /* 1 + 2^(1 - FLINT_BITS n) */
            NFLOAT_EXP(res) = 1;
            _nfloat_zero_limbs(NFLOAT_D(res), n);
            NFLOAT_D(res)[0] = 1;
            NFLOAT_D(res)[n - 1] |= UWORD(1) << (FLINT_BITS - 1);
            return GR_SUCCESS;
        }

        if (!up && dsgn)
        {
            /* 1 - 2^(-FLINT_BITS n) */
            slong i;
            NFLOAT_EXP(res) = 0;
            for (i = 0; i < n; i++)
                NFLOAT_D(res)[i] = UWORD_MAX;
            return GR_SUCCESS;
        }
    }

    NFLOAT_EXP(res) = 1;
    _nfloat_zero_limbs(NFLOAT_D(res), n - 1);
    NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
    return GR_SUCCESS;
}

int
_nfloat_set_tiny_rel_err(nfloat_ptr res, nfloat_srcptr x, int side, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx);
    int up;

    if (res != x)
        nfloat_set(res, x, ctx);

    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        return GR_SUCCESS;

    up = nfloat_should_round_up(NFLOAT_SGNBIT(x), ctx);

    if (up && side > 0)
    {
        /* the magnitude one ulp up */
        if (mpn_add_1(NFLOAT_D(res), NFLOAT_D(res), n, 1))
        {
            NFLOAT_D(res)[n - 1] = UWORD(1) << (FLINT_BITS - 1);
            NFLOAT_EXP(res)++;
            NFLOAT_HANDLE_OVERFLOW(res, ctx);
        }
    }
    else if (!up && side < 0)
    {
        /* the magnitude one ulp down */
        mpn_sub_1(NFLOAT_D(res), NFLOAT_D(res), n, 1);
        if (!LIMB_MSB_IS_SET(NFLOAT_D(res)[n - 1]))
        {
            /* was a power of two: all ones one binade lower */
            slong i;
            for (i = 0; i < n; i++)
                NFLOAT_D(res)[i] = UWORD_MAX;
            NFLOAT_EXP(res)--;
            NFLOAT_HANDLE_UNDERFLOW(res, ctx);
        }
    }

    return GR_SUCCESS;
}

ulong
_nfloat_elem_normalize(nn_ptr t, slong * te, nn_srcptr y, slong len, slong ue, ulong err, slong wn)
{
    slong sh, i;
    unsigned int c;

    while (y[len - 1] == 0)
        len--;

    c = flint_clz(y[len - 1]);
    sh = FLINT_BITS * len - c - FLINT_BITS * wn;    /* right shift */
    *te = ue + sh;

    if (sh >= 0)
    {
        slong ls = sh / FLINT_BITS;
        unsigned int bs = sh % FLINT_BITS;
        int dropped = 0;

        for (i = 0; i < ls; i++)
            dropped |= (y[i] != 0);

        if (bs == 0)
            _nfloat_copy_limbs(t, y + ls, wn);
        else
        {
            dropped |= ((y[ls] << (FLINT_BITS - bs)) != 0);
            for (i = 0; i < wn; i++)
                t[i] = (y[ls + i] >> bs) | ((ls + i + 1 < len) ? (y[ls + i + 1] << (FLINT_BITS - bs)) : 0);
        }

        if (sh >= FLINT_BITS)
            err = 1;
        else
            err = (err >> sh) + (sh != 0);

        return err + dropped;
    }
    else
    {
        slong ls = (-sh) / FLINT_BITS;
        unsigned int bs = (-sh) % FLINT_BITS;

        _nfloat_zero_limbs(t, ls);
        if (bs == 0)
            _nfloat_copy_limbs(t + ls, y, len);
        else
        {
            t[ls] = y[0] << bs;
            for (i = 1; i < len; i++)
                t[ls + i] = (y[i] << bs) | (y[i - 1] >> (FLINT_BITS - bs));
        }

        /* the error scales up by 2^-sh (UWORD_MAX signals overflow) */
        if (-sh >= FLINT_BITS - 2 || err > (UWORD(1) << (FLINT_BITS - 3 + sh)))
            return UWORD_MAX;
        return err << (-sh);
    }
}

ulong
_nfloat_elem_div(nn_ptr q, nn_srcptr a, ulong ea, nn_srcptr d, ulong ed, slong w)
{
    ulong t[2 * NFLOAT_ELEM_MAX_LIMBS];

    FLINT_ASSERT(LIMB_MSB_IS_SET(a[w - 1]));
    FLINT_ASSERT(LIMB_MSB_IS_SET(d[w - 1]));

    _nfloat_zero_limbs(t, w);
    _nfloat_copy_limbs(t + w, a, w);
    flint_mpn_divapprox(q, t, 2 * w, d, w);

    /* |A B^w / D - A' B^w / D'| <= B^w (ea / D + (A' / D') ed / D) with
       D >= B^w / 2 and A' / D' <= 2 (1 + 2^-20), plus one for the
       approximate quotient */
    return 2 * ea + 5 * ed + 2;
}

/* Approximate reals with relative error bounds ****************************/

#define RF_ERR_MAX (UWORD(1) << (FLINT_BITS - 4))

FLINT_FORCE_INLINE ulong
_rf_err_add(ulong a, ulong b)
{
    if (a >= RF_ERR_MAX || b >= RF_ERR_MAX)
        return UWORD_MAX;
    return a + b;
}

void
_nfloat_rf_set_mpn(_nfloat_rf_t r, nn_srcptr y, slong len, slong ue, ulong err, slong w)
{
    slong te;
    r->err = _nfloat_elem_normalize(r->d, &te, y, len, ue, err, w);
    r->exp = te + FLINT_BITS * w;
}

void
_nfloat_rf_mul(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 1];
    ulong err;
    slong exp;

    /* (A + ea)(B + eb) / B^w - AB / B^w <= ea + eb + ea eb / B^w, plus two
       for the high product */
    flint_mpn_mulhigh_n(t, a->d, b->d, w);
    err = _rf_err_add(_rf_err_add(a->err, b->err), 3);
    exp = a->exp + b->exp;

    if (!LIMB_MSB_IS_SET(t[w - 1]))
    {
        mpn_lshift(t, t, w, 1);
        err = _rf_err_add(err, err);
        exp--;
    }

    _nfloat_copy_limbs(r->d, t, w);
    r->exp = exp;
    r->err = err;
}

void
_nfloat_rf_div(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w)
{
    ulong q[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong exp;

    if (a->err >= RF_ERR_MAX || b->err >= RF_ERR_MAX)
    {
        r->err = UWORD_MAX;
        return;
    }

    err = _nfloat_elem_div(q, a->d, a->err, b->d, b->err, w);
    exp = a->exp - b->exp;

    /* Q = A B^w / D in (B^w / 2, 2 B^w) */
    if (q[w] != 0)
    {
        ulong lost = q[0] & 1;
        mpn_rshift(q, q, w + 1, 1);
        q[w - 1] |= (UWORD(1) << (FLINT_BITS - 1));
        err = (err >> 1) + 1 + lost;
        exp++;
    }

    _nfloat_copy_limbs(r->d, q, w);
    r->exp = exp;
    r->err = err;
}

void
_nfloat_rf_inv(_nfloat_rf_t r, const _nfloat_rf_t a, slong w)
{
    _nfloat_rf_t one;
    _nfloat_zero_limbs(one->d, w);
    one->d[w - 1] = UWORD(1) << (FLINT_BITS - 1);
    one->exp = 1;
    one->err = 0;
    _nfloat_rf_div(r, one, a, w);
}

/* T = (a, w) + or - (b, w) 2^-sh as w + 1 limbs, the shifted-out bits
   lost (one more unit) */
static ulong
_rf_aligned(nn_ptr T, nn_srcptr a, nn_srcptr b, slong sh, slong w, int sub, ulong eb)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 1];
    ulong err;
    slong ls = sh / FLINT_BITS, i;
    unsigned int bs = sh % FLINT_BITS;
    int dropped = 0;

    if (ls >= w)
    {
        _nfloat_zero_limbs(t, w);
        dropped = 1;
        err = 1;
    }
    else
    {
        for (i = 0; i < ls; i++)
            dropped |= (b[i] != 0);
        if (bs != 0)
        {
            dropped |= ((b[ls] << (FLINT_BITS - bs)) != 0);
            mpn_rshift(t, b + ls, w - ls, bs);
        }
        else
            _nfloat_copy_limbs(t, b + ls, w - ls);
        _nfloat_zero_limbs(t + w - ls, ls);
        err = (sh >= FLINT_BITS) ? 1 : (eb >> sh) + (sh != 0);
    }

    if (sub)
        T[w] = -mpn_sub_n(T, a, t, w);
    else
        T[w] = mpn_add_n(T, a, t, w);

    return err + dropped;
}

void
_nfloat_rf_add(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w)
{
    ulong T[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong te;

    if (a->exp < b->exp)
    {
        _nfloat_rf_add(r, b, a, w);
        return;
    }

    if (a->err >= RF_ERR_MAX || b->err >= RF_ERR_MAX)
    {
        r->err = UWORD_MAX;
        return;
    }

    err = a->err + _rf_aligned(T, a->d, b->d, a->exp - b->exp, w, 0, b->err);
    r->err = _nfloat_elem_normalize(r->d, &te, T, w + 1, a->exp - FLINT_BITS * w, err, w);
    r->exp = te + FLINT_BITS * w;
}

int
_nfloat_rf_sub(_nfloat_rf_t r, const _nfloat_rf_t a, const _nfloat_rf_t b, slong w)
{
    ulong T[NFLOAT_ELEM_MAX_LIMBS + 2];
    ulong err;
    slong te;
    int neg = 0;

    if (a->err >= RF_ERR_MAX || b->err >= RF_ERR_MAX)
    {
        r->err = UWORD_MAX;
        return 0;
    }

    if (a->exp < b->exp || (a->exp == b->exp && mpn_cmp(a->d, b->d, w) < 0))
    {
        const _nfloat_rf_struct * t = a;
        a = b;
        b = t;
        neg = 1;
    }

    err = a->err + _rf_aligned(T, a->d, b->d, a->exp - b->exp, w, 1, b->err);

    if (flint_mpn_zero_p(T, w + 1))
    {
        r->err = UWORD_MAX;
        return neg;
    }

    r->err = _nfloat_elem_normalize(r->d, &te, T, w + 1, a->exp - FLINT_BITS * w, err, w);
    r->exp = te + FLINT_BITS * w;
    return neg;
}

int
_nfloat_set_rf(nfloat_ptr res, const _nfloat_rf_t r, int sgnbit, slong w, gr_ctx_t ctx)
{
    ulong t[NFLOAT_ELEM_MAX_LIMBS + 2];

    if (r->err >= RF_ERR_MAX)
        return GR_UNABLE;

    _nfloat_copy_limbs(t, r->d, w);
    return _nfloat_set_mpn_err(res, t, w, r->exp - FLINT_BITS * w, r->err, r->err, sgnbit, ctx);
}

/* Ball evaluation ***********************************************************/

void
_nfloat_mp_real_view(mp_real_t v, nn_ptr buf, nfloat_srcptr x, gr_ctx_t ctx)
{
    slong n = NFLOAT_CTX_NLIMBS(ctx), e = NFLOAT_EXP(x), q, size;
    unsigned int r;
    nn_ptr d = buf;

    FLINT_ASSERT(!NFLOAT_IS_SPECIAL(x));

    /* |x| = M 2^(e - FLINT_BITS n) = (M 2^r) B^(q + 1 - (n + 1)) with
       e = FLINT_BITS q + r, 0 <= r < FLINT_BITS */
    q = (e >= 0) ? e / FLINT_BITS : -((-e + FLINT_BITS - 1) / FLINT_BITS);
    r = (unsigned int) (e - FLINT_BITS * q);

    if (r == 0)
    {
        _nfloat_copy_limbs(d, NFLOAT_D(x), n);
        size = n;
        v->exp = q;
    }
    else
    {
        d[n] = mpn_lshift(d, NFLOAT_D(x), n, r);
        size = n + 1;
        v->exp = q + 1;
    }

    /* exact values carry no low zero limbs */
    while (d[0] == 0)
    {
        d++;
        size--;
    }

    v->d = d;
    v->alloc = size;
    v->size = size;
    v->negative = NFLOAT_SGNBIT(x);
    v->err = 0;
}

int
_nfloat_set_mp_real(nfloat_ptr res, const mp_real_t x, slong shift, int * status, gr_ctx_t ctx)
{
    slong size = x->size, mbits;
    nn_ptr t;
    TMP_INIT;

    if (size == 0)
        return 0;

    /* the relative radius err / mid below 2^-(func_prec + guard bits) */
    if (x->err != 0)
    {
        mbits = FLINT_BITS * (size - 1) + FLINT_BIT_COUNT(x->d[size - 1]) - 1;
        if (mbits - (slong) FLINT_BIT_COUNT(x->err) < NFLOAT_CTX_FUNC_PREC(ctx) + NFLOAT_ELEM_GUARD_BITS)
            return 0;
    }

    TMP_START;
    t = TMP_ALLOC((size + 1) * sizeof(ulong));
    flint_mpn_copyi(t, x->d, size);
    *status = _nfloat_set_mpn_err(res, t, size, FLINT_BITS * (x->exp - size) + shift,
        x->err, x->err, x->negative, ctx);
    TMP_END;
    return 1;
}

int
_nfloat_eval_ball(nfloat_ptr res, _nfloat_ball_func_t f, nfloat_srcptr x, nfloat_srcptr y, int op, slong min_prec, gr_ctx_t ctx)
{
    ulong xbuf[NFLOAT_MAX_LIMBS + 1], ybuf[NFLOAT_MAX_LIMBS + 1];
    mp_real_t X, Y, R;
    slong prec, max_prec, shift;
    int status = GR_UNABLE, ok;

    _nfloat_mp_real_view(X, xbuf, x, ctx);
    if (y != NULL)
        _nfloat_mp_real_view(Y, ybuf, y, ctx);
    mp_real_init(R);

    max_prec = 16 * NFLOAT_CTX_PREC(ctx) + 100000;
    if (min_prec <= NFLOAT_TRIG_BALL_MAX_EXP / 4 + 64)
        max_prec = FLINT_MAX(max_prec, 2 * min_prec);
    prec = FLINT_MAX(min_prec, NFLOAT_CTX_FUNC_PREC(ctx) + NFLOAT_ELEM_GUARD_BITS + 32);

    for ( ; prec <= max_prec; prec *= 2)
    {
        shift = 0;
        ok = f(R, &shift, X, (y != NULL) ? Y : NULL, op, prec);

        if (ok < 0)
        {
            status = nfloat_nan(res, ctx);
            break;
        }

        if (ok == 0)
            continue;

        if (mp_real_is_zero(R))
        {
            status = nfloat_zero(res, ctx);
            break;
        }

        if (_nfloat_set_mp_real(res, R, shift, &status, ctx))
            break;
    }

    mp_real_clear(R);
    return status;
}

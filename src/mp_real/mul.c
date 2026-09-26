/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Products: mp_real_mul, mp_real_mul_ui, mp_real_mul_complex and the
   windowed residual _mp_real_submul_bounded. */

/* the windowed middle product beats the full product from this many
   kept limbs on (measured on x86-64: at 32 limbs the window costs 1.3
   times the full product, at 64 the two are even, at 256 the window
   is 30% faster) */
#ifndef MP_REAL_MULMID_CUTOFF
#define MP_REAL_MULMID_CUTOFF 64
#endif

/* the radius of a product: the cross terms |a| errb + |b| erra +
   erra errb, valid also for zero mantissas, magnitudes through the
   leading bits (see _mp_real_lead: bounding through the top limb alone
   overstates by a factor up to 2 for a top limb of 1, and that
   pessimism compounds across a chain of multiplications) */
static mp_real_bnd_t
_mp_real_mul_bnd(const mp_real_t a, const mp_real_t b)
{
    mp_real_bnd_t e = _mp_real_bnd_zero;

    if (b->err != 0)
    {
        if (a->size == 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, b->err, a->exp + b->exp - b->size));
        else
            e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(b->err, a,
                    (a->exp - 1) + b->exp - b->size));
    }
    if (a->err != 0)
    {
        if (b->size == 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, a->err, b->exp + a->exp - a->size));
        else
            e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(a->err, b,
                    (b->exp - 1) + a->exp - a->size));
        if (b->err != 0)
        {
            mp_real_bnd_t t;
            umul_ppmm(t.hi, t.lo, a->err, b->err);
            t.a = a->exp - a->size + b->exp - b->size;
            e = _mp_real_bnd_add(e, t);
        }
    }
    return e;
}

/* res = a * b truncated to about n limbs.  Full products are used
   when at most a few limbs longer than the kept part (the low limbs
   then round down into <= 1 ulp) and below MP_REAL_MULMID_CUTOFF kept
   limbs, otherwise flint_mpn_mulmid computes just the kept window
   plus one guard limb below, whose one-sided deficit of up to
   min(an, bn) + 2 ulps one limb above the window bottom is absorbed
   into the radius (and immediately truncated away by the error
   normalization). */
void
mp_real_mul(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
{
    slong sa = a->size, sb = b->size, ps, keep;
    slong rexp;
    int rneg, exact = (a->err == 0 && b->err == 0);
    mp_real_bnd_t e;
    nn_srcptr ad, bd;
    TMP_INIT;

    e = _mp_real_mul_bnd(a, b);

    if (sa == 0 || sb == 0)
    {
        _mp_real_zero_bnd(res, e);
        return;
    }

    ps = sa + sb;

    /* short exact operands whose product is kept whole: no windows,
       scratch or bound arithmetic (the leaves of binary splitting) */
    if (exact && ps <= 4 && ps <= n)
    {
        ulong t[4];
        slong lo, sz;

        if (ps == 2)
            umul_ppmm(t[1], t[0], a->d[0], b->d[0]);
        else if (sa >= sb)
            flint_mpn_mul(t, a->d, sa, b->d, sb);
        else
            flint_mpn_mul(t, b->d, sb, a->d, sa);

        sz = ps - (t[ps - 1] == 0);
        for (lo = 0; t[lo] == 0; lo++)
            ;
        mp_real_fit_length(res, sz - lo);
        res->d[0] = t[lo];
        if (sz - lo > 1) res->d[1] = t[lo + 1];
        if (sz - lo > 2) res->d[2] = t[lo + 2];
        if (sz - lo > 3) res->d[3] = t[lo + 3];
        res->size = sz - lo;
        res->exp = a->exp + b->exp - (ps - sz);
        res->negative = a->negative ^ b->negative;
        res->err = 0;
        return;
    }

    keep = FLINT_MIN(n, ps);
    if (!exact)
    {
        keep = FLINT_MIN(keep, _mp_real_acc(a, n));
        keep = FLINT_MIN(keep, _mp_real_acc(b, n));
    }
    keep = FLINT_MAX(keep, 1);

    rexp = a->exp + b->exp;
    rneg = a->negative ^ b->negative;

    /* an operand longer than keep + 2 limbs enters the product by its
       top keep + 2 limbs only: the dropped tail, below one unit of
       the lowest kept limb, moves the product by less than the other
       operand's magnitude (through its top limb, as above) times
       that unit */
    ad = a->d;
    bd = b->d;
    if (sa > keep + 2)
    {
        e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(1, b,
                (a->exp - (keep + 2)) + (b->exp - 1)));
        ad += sa - (keep + 2);
        sa = keep + 2;
    }
    if (sb > keep + 2)
    {
        e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(1, a,
                (b->exp - (keep + 2)) + (a->exp - 1)));
        bd += sb - (keep + 2);
        sb = keep + 2;
    }
    ps = sa + sb;

    if (ps <= keep + 5 || ps <= keep + keep / 16
        || keep < MP_REAL_MULMID_CUTOFF || (ad == bd && sa == sb))
    {
        /* full product (a square always: flint_mpn_sqr at two thirds
           of a product beats the window of a middle product, which
           costs 0.9 of one whatever its size), written straight into the destination
           buffer unless it would overwrite an operand; when its top
           limb is zero one more limb is kept, or the window would
           hold only keep - 1 significant ones */
        slong drop;
        int aliased = (res == a || res == b);
        nn_ptr t;
        TMP_START;
        if (aliased)
            t = TMP_ALLOC(ps * sizeof(ulong));
        else
        {
            mp_real_fit_length(res, ps);
            t = res->d;
        }

        if (ad == bd && sa == sb)
            flint_mpn_sqr(t, ad, sa);
        else if (sa >= sb)
            flint_mpn_mul(t, ad, sa, bd, sb);
        else
            flint_mpn_mul(t, bd, sb, ad, sa);

        keep = FLINT_MIN(keep + (t[ps - 1] == 0), ps);
        drop = ps - keep;

        if (drop > 0)
        {
            /* keep exactness when the dropped tail is zero */
            slong i;
            int inexact = 0;
            for (i = 0; i < drop && !inexact; i++)
                inexact = (t[i] != 0);
            if (inexact)
                e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, rexp - ps + drop));
        }

        if (aliased)
        {
            mp_real_fit_length(res, keep);
            flint_mpn_copyi(res->d, t + drop, keep);
        }
        else if (drop > 0)
            flint_mpn_copyi(res->d, res->d + drop, keep);
        res->size = keep;
        res->exp = rexp;
        res->negative = rneg;
        res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
        _mp_real_norm(res);
        _mp_real_apply_bnd(res, e);
        TMP_END;
    }
    else
    {
        /* windowed middle product: window [ps - keep - 2, ps), one
           limb more when the top limb of the product is zero (short of
           a carry from below).  The two guard limbs hold the deficit of
           the window, up to (min(sa, sb) + 2) B units of its bottom
           limb, below one ulp of the result (with one guard limb it was
           that many ulps: 1500 for operands of 1500 limbs) */
        slong zlo, zn;
        int aliased = (res == a || res == b);
        nn_ptr t;
        TMP_START;

        {
            ulong hi, lo;
            umul_ppmm(hi, lo, ad[sa - 1], bd[sb - 1]);
            (void) lo;
            keep = FLINT_MIN(keep + (hi == 0), ps);
        }
        zlo = FLINT_MAX(ps - keep - 2, 0);
        zn = ps - zlo;

        if (aliased)
            t = TMP_ALLOC(zn * sizeof(ulong));
        else
        {
            mp_real_fit_length(res, zn);
            t = res->d;
        }

        flint_mpn_mulmid(t, ad, sa, bd, sb, zlo, ps);

        /* deficit (below the true window) plus the discarded low
           limbs: < (min + 2) B + 1 ulps of the window bottom */
        e = _mp_real_bnd_add(e, _mp_real_bnd(0, FLINT_MIN(FLINT_MIN(sa, sb), zlo) + 2,
            rexp - ps + zlo + 1));
        e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, rexp - ps + zlo));

        if (aliased)
        {
            mp_real_fit_length(res, zn);
            flint_mpn_copyi(res->d, t, zn);
        }
        res->size = zn;
        res->exp = rexp;
        res->negative = rneg;
        res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
        _mp_real_norm(res);
        _mp_real_apply_bnd(res, e);
        TMP_END;
    }
}

/* The window of one part of a complex operand at the common bottom
   limb bot of the operand: the part's limbs above bot (a pointer into
   its mantissa when its bottom lies at or below bot, else the mantissa
   over zero padding written into buf), at least one limb */
static nn_srcptr
_mp_real_cwin(const mp_real_t x, slong bot, nn_ptr buf, slong * len)
{
    slong xb = x->exp - x->size, k;

    if (x->size == 0 || x->exp <= bot)
    {
        buf[0] = 0;
        *len = 1;
        return buf;
    }
    if (xb >= bot)
    {
        k = xb - bot;
        flint_mpn_zero(buf, k);
        flint_mpn_copyi(buf + k, x->d, x->size);
        *len = k + x->size;
        return buf;
    }
    *len = x->exp - bot;
    return x->d + (bot - xb);
}

/* install the signed exact product (z, |zl|) B^bot on res, dropping
   its limbs below the position cut (a value entirely below cut
   becomes a zero ball covering it), with the bound e */
static void
_mp_real_set_product(mp_real_t res, nn_srcptr z, slong zl, slong bot,
    slong cut, mp_real_bnd_t e)
{
    slong len = FLINT_ABS(zl), drop;

    while (len > 0 && z[len - 1] == 0)
        len--;
    drop = FLINT_MAX(0, cut - bot);
    if (drop >= len)
    {
        if (len > 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, bot + len));
        if (_mp_real_bnd_is_zero(e))
        {
            res->size = 0;
            res->negative = 0;
            res->exp = cut;
            res->err = 0;
        }
        else
            _mp_real_zero_bnd(res, e);
        return;
    }
    if (drop > 0)
    {
        slong i;
        int inexact = 0;
        for (i = 0; i < drop && !inexact; i++)
            inexact = (z[i] != 0);
        if (inexact)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, bot + drop));
    }
    mp_real_fit_length(res, len - drop);
    flint_mpn_copyi(res->d, z + drop, len - drop);
    res->size = len - drop;
    res->exp = bot + len;
    res->negative = (zl < 0);
    res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
    _mp_real_norm(res);
    _mp_real_apply_bnd(res, e);
}

/* |x| < B^mexp(x), for a zero mantissa through the radius */
static slong
_mp_real_mexp(const mp_real_t x)
{
    return x->size ? x->exp : x->exp + 1;
}

/* rr + i ri = (ar + i ai) (br + i bi) to about n limbs IN ONE FRAME:
   the precision is relative to the larger part of the result, both
   parts kept down to the same limb (a much smaller part keeps
   correspondingly fewer limbs, or becomes a zero ball inside its
   radius) -- the precision of a complex number is that of its
   modulus, and an accumulation like prod exp(i x_k) needs its small
   imaginary part to the absolute precision of the frame, not to its
   own relative precision (which would cost windows longer by the
   ratio of the parts).  Likewise each operand's parts are windowed
   to a common bottom keep + 2 limbs below the larger part.  One
   complex product of flint_mpn_mul_complex over the four windows
   (above the middle-product cutoff its high half only, by
   flint_mpn_mulhigh_n_complex: four forward and two inverse
   transforms in place of the twelve of four products).  The outputs
   may alias the inputs. */
void
mp_real_mul_complex(mp_real_t rr, mp_real_t ri, const mp_real_t ar,
    const mp_real_t ai, const mp_real_t br, const mp_real_t bi, slong n)
{
    slong topa, topb, T, fl, keep, bota, botb, cut, la[2], lb[2], zl;
    slong ea_r, ea_i, eb_r, eb_i, zr_len, zi_len;
    nn_srcptr wa[2], wb[2];
    nn_ptr buf, zr, zi;
    mp_real_bnd_t er, ei;
    TMP_INIT;

    ea_r = _mp_real_mexp(ar); ea_i = _mp_real_mexp(ai);
    eb_r = _mp_real_mexp(br); eb_i = _mp_real_mexp(bi);

    if ((ar->size == 0 && ai->size == 0) || (br->size == 0 && bi->size == 0))
    {
        /* a zero operand (possibly with a radius): the plain formulas */
        mp_real_t t1, t2;
        mp_real_init(t1);
        mp_real_init(t2);
        mp_real_mul(t1, ar, br, n);
        mp_real_mul(t2, ai, bi, n);
        mp_real_sub(t1, t1, t2, n);
        mp_real_mul(t2, ar, bi, n);
        mp_real_mul(ri, ai, br, n);
        mp_real_add(ri, ri, t2, n);
        mp_real_swap(rr, t1);
        mp_real_clear(t1);
        mp_real_clear(t2);
        return;
    }

    /* the frames: the tops of the operands and of the result, the
       noise floor of the result (the term p q with p inexact is noisy
       from mexp(q) + floor(p) up), and the kept limbs */
    topa = FLINT_MAX(ar->size ? ar->exp : WORD_MIN, ai->size ? ai->exp : WORD_MIN);
    topb = FLINT_MAX(br->size ? br->exp : WORD_MIN, bi->size ? bi->exp : WORD_MIN);
    T = topa + topb;
    fl = WORD_MIN;
    if (ar->err) fl = FLINT_MAX(fl, FLINT_MAX(eb_r, eb_i) + _mp_real_noise_floor(ar));
    if (ai->err) fl = FLINT_MAX(fl, FLINT_MAX(eb_r, eb_i) + _mp_real_noise_floor(ai));
    if (br->err) fl = FLINT_MAX(fl, FLINT_MAX(ea_r, ea_i) + _mp_real_noise_floor(br));
    if (bi->err) fl = FLINT_MAX(fl, FLINT_MAX(ea_r, ea_i) + _mp_real_noise_floor(bi));
    keep = n;
    if (fl != WORD_MIN)
        keep = FLINT_MIN(keep, T - fl + MP_REAL_ACC_GUARD);
    keep = FLINT_MAX(keep, 1);
    cut = T - keep;

    /* the windows: each operand's parts down to keep + 2 limbs below
       its top, never below the data */
    bota = topa - (keep + 2);
    botb = topb - (keep + 2);
    bota = FLINT_MAX(bota, FLINT_MIN(ar->size ? ar->exp - ar->size : WORD_MAX,
                                     ai->size ? ai->exp - ai->size : WORD_MAX));
    botb = FLINT_MAX(botb, FLINT_MIN(br->size ? br->exp - br->size : WORD_MAX,
                                     bi->size ? bi->exp - bi->size : WORD_MAX));
    la[0] = ar->size ? FLINT_MAX(ar->exp - bota, 1) : 1;
    la[1] = ai->size ? FLINT_MAX(ai->exp - bota, 1) : 1;
    lb[0] = br->size ? FLINT_MAX(br->exp - botb, 1) : 1;
    lb[1] = bi->size ? FLINT_MAX(bi->exp - botb, 1) : 1;

    /* the bounds: the cross terms of the four products, and the
       truncation of a part below its window (under one unit there)
       times the magnitude of the parts it multiplies */
    er = _mp_real_bnd_add(_mp_real_mul_bnd(ar, br), _mp_real_mul_bnd(ai, bi));
    ei = _mp_real_bnd_add(_mp_real_mul_bnd(ar, bi), _mp_real_mul_bnd(ai, br));
    if (ar->size && ar->exp - ar->size < bota)
    {
        if (br->size) er = _mp_real_bnd_add(er, _mp_real_bnd_mul_mag(1, br, bota + br->exp - 1));
        if (bi->size) ei = _mp_real_bnd_add(ei, _mp_real_bnd_mul_mag(1, bi, bota + bi->exp - 1));
    }
    if (ai->size && ai->exp - ai->size < bota)
    {
        if (bi->size) er = _mp_real_bnd_add(er, _mp_real_bnd_mul_mag(1, bi, bota + bi->exp - 1));
        if (br->size) ei = _mp_real_bnd_add(ei, _mp_real_bnd_mul_mag(1, br, bota + br->exp - 1));
    }
    if (br->size && br->exp - br->size < botb)
    {
        if (ar->size) er = _mp_real_bnd_add(er, _mp_real_bnd_mul_mag(1, ar, botb + ar->exp - 1));
        if (ai->size) ei = _mp_real_bnd_add(ei, _mp_real_bnd_mul_mag(1, ai, botb + ai->exp - 1));
    }
    if (bi->size && bi->exp - bi->size < botb)
    {
        if (ai->size) er = _mp_real_bnd_add(er, _mp_real_bnd_mul_mag(1, ai, botb + ai->exp - 1));
        if (ar->size) ei = _mp_real_bnd_add(ei, _mp_real_bnd_mul_mag(1, ar, botb + ar->exp - 1));
    }

    TMP_START;
    if (keep >= MP_REAL_MULMID_CUTOFF
        && 4 * FLINT_MIN(FLINT_MAX(la[0], la[1]), FLINT_MAX(lb[0], lb[1]))
            >= 3 * FLINT_MAX(FLINT_MAX(la[0], la[1]), FLINT_MAX(lb[0], lb[1])))
    {
        /* the high half only, for balanced windows (a short exact
           operand is cheaper multiplied in full): the four parts
           brought to a common length nn by zero padding below
           (shifting an operand's scale down by its padding, and
           padded further when the window would not reach the cut),
           the limbs [nn, 2 nn] of the exact product come back within
           4 units of their lowest limb, and the discarded low half
           is under one more */
        slong nn = FLINT_MAX(FLINT_MAX(la[0], la[1]), FLINT_MAX(lb[0], lb[1])) + 2;
        slong ka, kb, scale;
        int sr, si;
        nn_ptr pa[2], pb[2];

        scale = (bota + botb) - 2 * nn + FLINT_MAX(la[0], la[1]) + FLINT_MAX(lb[0], lb[1]) + nn;
        if (scale > cut)
            nn += scale - cut;
        ka = nn - FLINT_MAX(la[0], la[1]);
        kb = nn - FLINT_MAX(lb[0], lb[1]);
        scale = (bota - ka) + (botb - kb) + nn;

        buf = TMP_ALLOC((4 * nn + 2 * (nn + 2)) * sizeof(ulong));
        pa[0] = buf; pa[1] = buf + nn; pb[0] = buf + 2 * nn; pb[1] = buf + 3 * nn;
        zr = buf + 4 * nn;
        zi = zr + nn + 2;
        {
            slong l;
            nn_srcptr w;
            w = _mp_real_cwin(ar, bota, zr, &l);
            flint_mpn_zero(pa[0], ka); flint_mpn_copyi(pa[0] + ka, w, l); flint_mpn_zero(pa[0] + ka + l, nn - ka - l);
            w = _mp_real_cwin(ai, bota, zr, &l);
            flint_mpn_zero(pa[1], ka); flint_mpn_copyi(pa[1] + ka, w, l); flint_mpn_zero(pa[1] + ka + l, nn - ka - l);
            w = _mp_real_cwin(br, botb, zr, &l);
            flint_mpn_zero(pb[0], kb); flint_mpn_copyi(pb[0] + kb, w, l); flint_mpn_zero(pb[0] + kb + l, nn - kb - l);
            w = _mp_real_cwin(bi, botb, zr, &l);
            flint_mpn_zero(pb[1], kb); flint_mpn_copyi(pb[1] + kb, w, l); flint_mpn_zero(pb[1] + kb + l, nn - kb - l);
        }

        flint_mpn_mulhigh_n_complex(zr, &sr, zi, &si,
            pa[0], ar->negative, pa[1], ai->negative,
            pb[0], br->negative, pb[1], bi->negative, nn);

        er = _mp_real_bnd_add(er, _mp_real_bnd(0, 5, scale));
        ei = _mp_real_bnd_add(ei, _mp_real_bnd(0, 5, scale));
        zl = nn + 1;
        _mp_real_set_product(rr, zr, sr ? -zl : zl, scale, cut, er);
        _mp_real_set_product(ri, zi, si ? -zl : zl, scale, cut, ei);
        TMP_END;
        return;
    }

    zl = FLINT_MAX(la[0], la[1]) + FLINT_MAX(lb[0], lb[1]) + 1;
    buf = TMP_ALLOC((la[0] + la[1] + lb[0] + lb[1] + 2 * zl) * sizeof(ulong));
    wa[0] = _mp_real_cwin(ar, bota, buf, &la[0]);
    wa[1] = _mp_real_cwin(ai, bota, buf + la[0], &la[1]);
    wb[0] = _mp_real_cwin(br, botb, buf + la[0] + la[1], &lb[0]);
    wb[1] = _mp_real_cwin(bi, botb, buf + la[0] + la[1] + lb[0], &lb[1]);
    zr = buf + la[0] + la[1] + lb[0] + lb[1];
    zi = zr + zl;

    flint_mpn_mul_complex(zr, &zr_len, zi, &zi_len,
        wa[0], la[0], ar->negative, wa[1], la[1], ai->negative,
        wb[0], lb[0], br->negative, wb[1], lb[1], bi->negative);

    _mp_real_set_product(rr, zr, zr_len, bota + botb, cut, er);
    _mp_real_set_product(ri, zi, zi_len, bota + botb, cut, ei);
    TMP_END;
}

void
mp_real_mul_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)
{
    mp_real_bnd_t e;
    nn_ptr t;
    ulong cy;
    slong sa = a->size, keep, drop;

    if (c == 0)
    {
        mp_real_zero(res);
        return;
    }

    e = _mp_real_bnd_zero;
    if (a->err != 0)
    {
        umul_ppmm(e.hi, e.lo, a->err, c);
        e.a = (sa == 0) ? a->exp : a->exp - a->size;
    }

    if (sa == 0)
    {
        _mp_real_zero_bnd(res, e);
        return;
    }

    /* works in the destination buffer, in place when aliased
       (mp_real_fit_length preserves contents across reallocation) */
    mp_real_fit_length(res, sa + 1);
    t = res->d;
    cy = mpn_mul_1(t, a->d, sa, c);
    t[sa] = cy;

    {
        slong anc = a->exp - a->size;  /* ulp of t */
        slong ts = sa + (cy != 0);     /* top-normalized length */

        keep = FLINT_MIN(n + 1, ts);
        drop = ts - keep;
        if (drop > 0)
        {
            slong i;
            int inexact = 0;
            for (i = 0; i < drop && !inexact; i++)
                inexact = (t[i] != 0);
            if (inexact)
                e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, anc + drop));
            flint_mpn_copyi(t, t + drop, keep);
        }

        res->size = keep;
        res->exp = anc + ts;
        res->negative = a->negative;
        res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
        _mp_real_norm(res);
        _mp_real_apply_bnd(res, e);
    }
}

/* res = a - b c to n limbs, GIVEN that |a - b c| < B^E (the caller's
   knowledge, typically from the analysis of an iteration: a residual
   1 - x z^k with z accurate to some precision).  Only the window of
   the product b c between the weights B^(E - n - 2) and B^E is
   computed (flint_mpn_mulmid, from the limbs of b and c that reach
   it), against which the same window of a is taken: the difference
   modulo B^(E+1) is the residual in two's complement, its limb of
   weight B^E being 0 or B - 1 by the assumption, which gives the
   sign.  The two limbs below the retained n absorb the deficit of the
   middle product and the dropped tails of the operands, so that the
   truncation costs one unit of the result's bottom limb; the radii of
   the operands enter as in mp_real_mul.  With the assumption violated
   the result is wrong (the window wraps), so it must be rigorous: a
   ball computed by an ordinary mp_real_mul and mp_real_sub would have been
   correct at any magnitude, at the cost of a high product of b c that
   is known to cancel. */
void
_mp_real_submul_bounded(mp_real_t res, const mp_real_t a, const mp_real_t b,
    const mp_real_t c, slong E, slong n)
{
    slong sb = b->size, sc = c->size, sa = a->size;
    slong wlo = E - n - 2, wn = n + 3, base, zlo, zhi, i, lo, hi;
    nn_srcptr bd, cd;
    nn_ptr t;
    mp_real_bnd_t e;
    int neg;
    TMP_INIT;

    /* the radii: err_a + |b| err_c + |c| err_b + err_b err_c, through
       the top limbs as in mp_real_mul */
    e = _mp_real_bnd_of(a);
    if (c->err != 0)
    {
        if (sb == 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, c->err, b->exp + c->exp - c->size));
        else
            e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(c->err, b,
                    (b->exp - 1) + c->exp - c->size));
    }
    if (b->err != 0)
    {
        if (sc == 0)
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, b->err, c->exp + b->exp - b->size));
        else
            e = _mp_real_bnd_add(e, _mp_real_bnd_mul_mag(b->err, c,
                    (c->exp - 1) + b->exp - b->size));
        if (c->err != 0)
        {
            mp_real_bnd_t t2;
            umul_ppmm(t2.hi, t2.lo, b->err, c->err);
            t2.a = b->exp - b->size + c->exp - c->size;
            e = _mp_real_bnd_add(e, t2);
        }
    }

    TMP_START;
    t = TMP_ALLOC(2 * wn * sizeof(ulong));

    /* the window of |b c|: limb i of the product has weight
       B^(base + i); the operands enter by the limbs that reach the
       window bottom (a limb of b of weight w meets all of c below
       B^(c.exp), so those with w + c.exp <= wlo - 1 are dropped, the
       tail being below B^(wlo - 1) times |c| < B^(wlo)... within the
       two guard limbs) */
    flint_mpn_zero(t, wn);
    if (sb > 0 && sc > 0)
    {
        slong bskip = FLINT_MAX(0, (wlo - 1) - c->exp - (b->exp - sb));
        slong cskip = FLINT_MAX(0, (wlo - 1) - b->exp - (c->exp - sc));
        slong sbn = sb - FLINT_MIN(bskip, sb);
        slong scn = sc - FLINT_MIN(cskip, sc);

        if (sbn > 0 && scn > 0)
        {
            bd = b->d + (sb - sbn);
            cd = c->d + (sc - scn);
            base = (b->exp - sbn) + (c->exp - scn);
            /* window [wlo, E] in product limb indices, clipped to the
               product's sbn + scn limbs */
            lo = wlo - base;
            hi = E + 1 - base;
            zlo = FLINT_MAX(lo, 0);
            zhi = FLINT_MIN(hi, sbn + scn);
            if (zlo < zhi)
            {
                if (bd == cd && sbn == scn)
                {
                    /* a square: in full by flint_mpn_sqr, cheaper than
                       any window of a middle product */
                    nn_ptr sq = TMP_ALLOC(2 * sbn * sizeof(ulong));
                    flint_mpn_sqr(sq, bd, sbn);
                    flint_mpn_copyi(t + wn + (zlo - lo), sq + zlo, zhi - zlo);
                }
                else if (sbn >= scn)
                    flint_mpn_mulmid(t + wn + (zlo - lo), bd, sbn, cd, scn,
                        zlo, zhi);
                else
                    flint_mpn_mulmid(t + wn + (zlo - lo), cd, scn, bd, sbn,
                        zlo, zhi);
                /* limbs of the window outside the product are zero */
                flint_mpn_zero(t + wn, zlo - lo);
                flint_mpn_zero(t + wn + (zhi - lo), wn - (zhi - lo));
                /* the deficit of the middle product (a lower
                   approximation) and the dropped tails: below
                   (min(sbn, scn) + 2) B + 2 units of B^wlo, absorbed by
                   the guard limbs into one unit of the result's ulp
                   together with the truncation below */
                (void) 0;
            }
            else
                flint_mpn_zero(t + wn, wn);
        }
        else
            flint_mpn_zero(t + wn, wn);
    }
    else
        flint_mpn_zero(t + wn, wn);

    /* the window of |a| */
    if (sa > 0)
    {
        slong abase = a->exp - sa;
        lo = wlo - abase;
        hi = E + 1 - abase;
        zlo = FLINT_MAX(lo, 0);
        zhi = FLINT_MIN(hi, sa);
        if (zlo < zhi)
            flint_mpn_copyi(t + (zlo - lo), a->d + zlo, zhi - zlo);
    }

    /* (-1)^a.neg |a| - (-1)^(b.neg + c.neg) |b c| modulo B^(E+1) */
    {
        int pneg = b->negative ^ c->negative;
        if (a->negative)
            mpn_neg(t, t, wn);
        if (pneg)
            mpn_add_n(t, t, t + wn, wn);
        else
            mpn_sub_n(t, t, t + wn, wn);
    }

    /* the sign from the control limb, |residual| < B^E */
    neg = (t[wn - 1] != 0);
    FLINT_ASSERT(t[wn - 1] == 0 || t[wn - 1] == UWORD_MAX);
    if (neg)
        mpn_neg(t, t, wn);

    /* the retained n limbs, of weights B^(wlo + 2) .. B^(E - 1); the
       guard limbs, deficits and tails are below 3 units of B^(wlo+2)
       (the deficit (min + 2) B + 2 units of B^wlo is below B^2 units
       for min + 2 < B) */
    mp_real_fit_length(res, n);
    flint_mpn_copyi(res->d, t + 2, n);
    for (i = n; i > 0 && res->d[i - 1] == 0; i--)
        ;
    res->size = i;
    res->exp = wlo + 2 + i;
    res->negative = neg;
    e = _mp_real_bnd_add(e, _mp_real_bnd(0, 3, wlo + 2));
    res->err = 1;
    if (res->size == 0)
        _mp_real_zero_bnd(res, e);
    else
        _mp_real_apply_bnd(res, e);

    TMP_END;
}

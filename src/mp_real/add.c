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

/* Sums and differences, and the word multiply-accumulates. */

/* zero a few limbs without a memset call */
FLINT_FORCE_INLINE void
_mp_real_zerow(nn_ptr p, slong k)
{
    if (k > 4)
        flint_mpn_zero(p, k);
    else
    {
        switch (k)
        {
            case 4: p[3] = 0; /* fall through */
            case 3: p[2] = 0; /* fall through */
            case 2: p[1] = 0; /* fall through */
            case 1: p[0] = 0; /* fall through */
            default: break;
        }
    }
}

FLINT_FORCE_INLINE ulong
_mp_real_addw(nn_ptr t, slong tn, nn_srcptr b, slong bn)
{
    slong i;
    ulong cy = 0, x, y;
    if (tn > 4)
        return mpn_add(t, t, tn, b, bn);
    for (i = 0; i < bn; i++)
    {
        x = t[i] + cy;
        cy = (x < cy);
        y = x + b[i];
        cy += (y < x);
        t[i] = y;
    }
    for (; i < tn && cy; i++)
    {
        t[i]++;
        cy = (t[i] == 0);
    }
    return cy;
}

FLINT_FORCE_INLINE ulong
_mp_real_subw(nn_ptr t, slong tn, nn_srcptr b, slong bn)
{
    slong i;
    ulong bw = 0, x, y;
    if (tn > 4)
        return mpn_sub(t, t, tn, b, bn);
    for (i = 0; i < bn; i++)
    {
        x = t[i];
        y = x - bw;
        bw = (y > x);
        x = y - b[i];
        bw += (x > y);
        t[i] = x;
    }
    for (; i < tn && bw; i++)
    {
        bw = (t[i] == 0);
        t[i]--;
    }
    return bw;
}


/* Exact short operands (at most two limbs each, exponents within two
   limbs, so that the exact sum has at most five limbs), formed in a
   fixed window of registers: no bounds, no truncation, no calls.
   Returns 0 when the operands do not qualify. */
FLINT_FORCE_INLINE int
_mp_real_add_short(mp_real_t res, const mp_real_t a, int aneg, const mp_real_t b,
    int bneg, slong n)
{
    ulong w[5], v[5];
    slong sa = a->size, sb = b->size, top, bot, oa, ob, i, k;
    int rneg;

    if ((a->err | b->err) != 0 || (ulong) (sa - 1) > 1 || (ulong) (sb - 1) > 1
        || (ulong) (a->exp - b->exp + 2) > 4 || n < 3)
        return 0;

    top = FLINT_MAX(a->exp, b->exp) + 1;
    bot = FLINT_MIN(a->exp - a->size, b->exp - b->size);
    oa = (a->exp - a->size) - bot;
    ob = (b->exp - b->size) - bot;
    /* top - bot <= 5 */

    w[0] = w[1] = w[2] = w[3] = w[4] = 0;
    v[0] = v[1] = v[2] = v[3] = v[4] = 0;
    w[oa] = a->d[0];
    if (sa == 2)
        w[oa + 1] = a->d[1];
    v[ob] = b->d[0];
    if (sb == 2)
        v[ob + 1] = b->d[1];

    rneg = aneg;
    if (aneg == bneg)
    {
        ulong cy = 0;
        for (i = 0; i < 5; i++)
        {
            ulong x = w[i] + cy;
            cy = (x < cy);
            x += v[i];
            cy += (x < v[i]);
            w[i] = x;
        }
    }
    else
    {
        ulong bw = 0;
        /* |w| >= |v|: the higher window wins, comparing from the top */
        int ge = 1;
        for (i = 4; i >= 0; i--)
        {
            if (w[i] != v[i])
            {
                ge = (w[i] > v[i]);
                break;
            }
        }
        if (!ge)
        {
            for (i = 0; i < 5; i++)
            {
                ulong x = w[i];
                w[i] = v[i];
                v[i] = x;
            }
            rneg = bneg;
        }
        for (i = 0; i < 5; i++)
        {
            ulong x = w[i] - bw;
            bw = (x > w[i]);
            bw += (x < v[i]);
            w[i] = x - v[i];
        }
    }

    /* strip zero limbs at both ends (exact) */
    k = top - bot;
    while (k > 0 && w[k - 1] == 0)
        k--;
    if (k == 0)
    {
        res->size = 0;
        res->negative = 0;
        res->exp = 0;
        res->err = 0;
        return 1;
    }
    for (i = 0; w[i] == 0; i++)
        ;
    mp_real_fit_length(res, k - i);
    res->size = k - i;
    res->exp = bot + k;
    res->negative = rneg;
    res->err = 0;
    switch (k - i)
    {
        case 5: res->d[4] = w[i + 4]; /* fall through */
        case 4: res->d[3] = w[i + 3]; /* fall through */
        case 3: res->d[2] = w[i + 2]; /* fall through */
        case 2: res->d[1] = w[i + 1]; /* fall through */
        default: res->d[0] = w[i];
    }
    return 1;
}

/* |a| >= |b| for nonzero mantissas (top-normalized: the exponents
   decide, then the aligned top limbs, then the longer operand) */
FLINT_FORCE_INLINE int
_mp_real_mag_ge(const mp_real_t a, const mp_real_t b)
{
    slong m;
    int c;
    if (a->exp != b->exp)
        return a->exp > b->exp;
    m = FLINT_MIN(a->size, b->size);
    c = mpn_cmp(a->d + a->size - m, b->d + b->size - m, m);
    if (c != 0)
        return c > 0;
    return a->size >= b->size;
}

/* signed addition core: res = (-1)^aneg |a| + (-1)^bneg |b|.  The
   window of the first operand is built in the destination buffer (in
   place when the result aliases it) and the second is added to or
   subtracted from it there, with no temporary: the sum being
   symmetric, the operands are swapped so that the result never
   aliases the second one, and a difference is taken from the larger
   operand (decided up front, so that no negation pass is needed)
   except when the result aliases the smaller, where the window is
   subtracted from the second operand in place. */
static void
_mp_real_add(mp_real_t res, const mp_real_t a, int aneg, const mp_real_t b,
    int bneg, slong n)
{
    slong top, bot, span, oa, ob;
    mp_real_bnd_t e;
    nn_ptr t;
    int rneg, rev = 0, trunc_a = 0, trunc_b = 0;

    if (_mp_real_add_short(res, a, aneg, b, bneg, n))
        return;

    if ((a->err | b->err) == 0)
        e = _mp_real_bnd_zero;
    else
        e = _mp_real_bnd_add(_mp_real_bnd_of(a), _mp_real_bnd_of(b));

    if (a->size == 0 || b->size == 0)
    {
        if (a->size == 0 && b->size == 0)
        {
            _mp_real_zero_bnd(res, e);
            return;
        }
        /* e already contains both radii; apply_bnd overwrites err */
        if (a->size == 0)
        {
            mp_real_set(res, b);
            res->negative = bneg;
        }
        else
        {
            mp_real_set(res, a);
            res->negative = aneg;
        }
        _mp_real_apply_bnd(res, e);
        return;
    }

    if (a == b)
    {
        if (aneg == bneg)
        {
            mp_real_mul_ui(res, a, 2, n);
            res->negative = aneg;
        }
        else
            _mp_real_zero_bnd(res, e);
        return;
    }

    if (res == b)
    {
        const mp_real_struct * s = a; a = b; b = s;
        rneg = aneg; aneg = bneg; bneg = rneg;
    }

    if (aneg != bneg && !_mp_real_mag_ge(a, b))
    {
        if (res != a)
        {
            const mp_real_struct * s = a; a = b; b = s;
            rneg = aneg; aneg = bneg; bneg = rneg;
        }
        else
            rev = 1;
    }

    top = FLINT_MAX(a->exp, b->exp) + 1;
    bot = FLINT_MIN(a->exp - a->size, b->exp - b->size);
    bot = FLINT_MAX(bot, top - (n + 2));
    /* nothing below an inexact operand's noise floor (see
       _mp_real_noise_floor): the other operand's tail there is dropped
       for one unit of the floor, and an in-place operand never moves */
    if (a->err != 0)
        bot = FLINT_MAX(bot, _mp_real_noise_floor(a));
    if (b->err != 0)
        bot = FLINT_MAX(bot, _mp_real_noise_floor(b));

    span = top - bot;
    oa = (a->exp - a->size) - bot;

    mp_real_fit_length(res, span);
    t = res->d;

    if (oa >= 0)
    {
        /* a fits fully: oa + a->size = a->exp - bot <= span - 1 */
        if (res == a)
        {
            if (oa > 0)
                memmove(t + oa, t, a->size * sizeof(ulong));
        }
        else
            flint_mpn_copyi(t + oa, a->d, a->size);
        if (oa > 0)
            _mp_real_zerow(t, oa);
        _mp_real_zerow(t + oa + a->size, span - oa - a->size);
    }
    else if (a->size + oa > 0)
    {
        /* low limbs of a truncated; copying downward is safe in
           place (ascending copy, dst below src) */
        flint_mpn_copyi(t, (res == a ? t : a->d) - oa, a->size + oa);
        _mp_real_zerow(t + a->size + oa, span - a->size - oa);
        trunc_a = 1;
    }
    else
    {
        trunc_a = 1;
        _mp_real_zerow(t, span);
    }

    rneg = aneg;

    /* combine |b| */
    ob = (b->exp - b->size) - bot;
    {
        nn_srcptr bd = b->d;
        slong bs = b->size;
        if (ob < 0)
        {
            bd -= ob;
            bs += ob;
            ob = 0;
            trunc_b = 1;
        }
        if (bs > 0)
        {
            if (aneg == bneg)
            {
                ulong cy = _mp_real_addw(t + ob, span - ob, bd, bs);
                FLINT_ASSERT(cy == 0);  /* top slot reserved */
                (void) cy;
            }
            else if (!rev)
            {
                ulong bw = _mp_real_subw(t + ob, span - ob, bd, bs);
                FLINT_ASSERT(bw == 0);  /* |a| >= |b| */
                (void) bw;
            }
            else
            {
                /* t = |b| window - t, |b| >= |a|: the limbs of t above
                   the window of b are zero, so the difference is
                   confined to the window and the limbs below it */
                ulong bw = 0;
                if (ob > 0)
                    bw = mpn_neg(t, t, ob);
                bw = mpn_sub_n(t + ob, bd, t + ob, bs)
                   + mpn_sub_1(t + ob, t + ob, bs, bw);
                FLINT_ASSERT(bw == 0);
                (void) bw;
                rneg = bneg;
            }
        }
        else
        {
            /* b entirely below the window (then not the larger) */
            FLINT_ASSERT(!rev);
            trunc_b = 1;
        }
    }

    if (trunc_a)
        e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, bot));
    if (trunc_b)
        e = _mp_real_bnd_add(e, _mp_real_bnd(0, 1, bot));

    /* the reserved carry slot is usually empty */
    if (t[span - 1] == 0)
    {
        span--;
        top--;
    }
    res->size = span;
    res->exp = top;
    res->negative = rneg;
    res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
    _mp_real_norm(res);
    if (res->size == 0)
        _mp_real_zero_bnd(res, e);
    else
        _mp_real_apply_bnd(res, e);
}

void
mp_real_add(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
{
    _mp_real_add(res, a, a->negative, b, b->negative, n);
}

void
mp_real_sub(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
{
    _mp_real_add(res, a, a->negative, b, !b->negative, n);
}

/* res = x + (-1)^sub y c.  In place (res == x) whenever the product
   lands inside x's window -- y c below x's top limb, x no longer than
   the precision asks, for a difference |y c| < B^(x.exp - 1) <= |x|
   so that the sign cannot flip, and y's bottom at or above x's, or
   anywhere below it when x is inexact (its noise floor is its bottom
   limb, see _mp_real_noise_floor) -- by mpn_addmul_1 resp.
   mpn_submul_1 over y's limbs from x's bottom up and a carry or
   borrow that stops when it runs out: O(|y|), whatever the size of x.
   The dropped tail of y, below one unit of x's bottom limb, is below
   c units there after the multiplication: the high word of c times
   its top limb enters as a carry, and the rest is under 2 units.  The
   radius is x's plus c times y's (exact products of words), added as
   in mp_real_add.  Otherwise (and when res aliases y) y c is formed by
   mp_real_mul_ui and added by mp_real_add. */
static void
_mp_real_addmul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c,
    slong n, int sub)
{
    slong xbot, ybot, off, xs, ys, drop;
    int ysgn;
    mp_real_bnd_t e;

    if (c == 0 || mp_real_is_zero(y))
    {
        mp_real_set(res, x);
        return;
    }

    ysgn = y->negative ^ sub;
    xs = x->size;
    ys = y->size;
    xbot = x->exp - xs;
    ybot = y->exp - ys;
    drop = xbot - ybot;         /* limbs of y below x's bottom */

    if (res != y && xs > 0 && ys > 0 && xs <= n + 2
        && (drop <= 0 || (x->err != 0 && drop < ys))
        && y->exp + 1 <= x->exp
        && (ysgn == x->negative || y->exp + 1 < x->exp))
    {
        nn_ptr d;
        nn_srcptr yd = y->d;
        ulong cy, cin = 0;

        if (res != x)
            mp_real_set(res, x);
        d = res->d;

        e = _mp_real_bnd_of(res);
        if (drop > 0)
        {
            /* c y_low = hi + (lo + c rest) / B with rest < 1 in units
               of B^(xbot - 1): [hi, hi + 2) units of B^xbot */
            ulong lo;
            umul_ppmm(cin, lo, yd[drop - 1], c);
            (void) lo;
            yd += drop;
            ys -= drop;
            off = 0;
            e = _mp_real_bnd_add(e, _mp_real_bnd(0, 2, xbot));
        }
        else
            off = -drop;

        if (ysgn == res->negative)
        {
            cy = mpn_addmul_1(d + off, yd, ys, c);
            if (cin != 0)
                cy += mpn_add_1(d + off, d + off, ys, cin);
            if (cy != 0)
                cy = mpn_add_1(d + off + ys, d + off + ys, xs - off - ys, cy);
            if (cy != 0)
            {
                mp_real_fit_length(res, xs + 1);
                res->d[xs] = cy;
                res->size = xs + 1;
                res->exp++;
            }
        }
        else
        {
            cy = mpn_submul_1(d + off, yd, ys, c);
            if (cin != 0)
                cy += mpn_sub_1(d + off, d + off, ys, cin);
            if (cy != 0)
                cy = mpn_sub_1(d + off + ys, d + off + ys, xs - off - ys, cy);
            FLINT_ASSERT(cy == 0);
        }

        /* c times y's radius at y's anchor */
        if (y->err != 0)
        {
            mp_real_bnd_t t;
            umul_ppmm(t.hi, t.lo, y->err, c);
            t.a = ybot;
            e = _mp_real_bnd_add(e, t);
        }
        res->err = !_mp_real_bnd_is_zero(e);    /* norm must know inexactness */
        _mp_real_norm(res);
        if (res->size == 0)
            _mp_real_zero_bnd(res, e);
        else
            _mp_real_apply_bnd(res, e);
        return;
    }

    {
        mp_real_t t;
        mp_real_init(t);
        mp_real_mul_ui(t, y, c, n);
        if (sub)
            mp_real_sub(res, x, t, n);
        else
            mp_real_add(res, x, t, n);
        mp_real_clear(t);
    }
}

void
mp_real_addmul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c,
    slong n)
{
    _mp_real_addmul_ui(res, x, y, c, n, 0);
}

void
mp_real_submul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c,
    slong n)
{
    _mp_real_addmul_ui(res, x, y, c, n, 1);
}

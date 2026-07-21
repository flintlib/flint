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
#include "mpn_extras.h"
#include "fixed.h"
#include "arb.h"

/* ANCHORED ERROR BOUNDS.  During the composition of an error bound
   the natural anchors of the contributions (the ulp positions of the
   operands, of the output, of a truncation) differ by amounts far
   exceeding the double exponent range, so bounds are carried as
   pairs

       e = v * B^a,    v a double >= 0, a a limb exponent (slong),

   with the invariant (established by fberr_reduce) that v == 0.0 or
   1.0 <= v < 2^128.  All combining operations round UP:

   - reduction v >= 2^128 -> v * 2^-64 + 1 per limb: the scaling is an
     exact power-of-two multiplication and the discarded fraction is
     covered by the added ulp;
   - rounding any nonzero bound up to at least one ulp of its anchor
     costs nothing real and makes additions safe: when adding a bound
     whose anchor sits d >= 1 limbs lower, the scaled term
     ldexp(v, -64 d) may round or (for d >= 17, where v * 2^-1088
     with v < 2^128 drops below 2^-960 and eventually to subnormals
     and 0) underflow entirely, but the final multiplication by
     FBALL_EPS adds at least 2^-24 * max(v) >= 2^-24 to the sum,
     dominating everything that was lost. */

/* powers of the limb radix as doubles (exact) */
#if FLINT_BITS == 64
#define FBALL_D_B 0x1p64
#define FBALL_D_BINV 0x1p-64
#define FBALL_D_SQRTB 0x1p32
#define FBALL_D_SQRTBINV 0x1p-32
#define FBALL_LGB 6
#else
#define FBALL_D_B 0x1p32
#define FBALL_D_BINV 0x1p-32
#define FBALL_D_SQRTB 0x1p16
#define FBALL_D_SQRTBINV 0x1p-16
#define FBALL_LGB 5
#endif

typedef struct { double v; slong a; } fberr_t;

static void
fberr_reduce(fberr_t * x)
{
    if (x->v == 0.0)
        return;
    while (x->v >= 0x1p128)
    {
        x->v = x->v * FBALL_D_BINV + 1.0;
        x->a++;
    }
    /* keep 1 <= v < 2^128 by REBASING the anchor rather than
       flooring: v B^a is preserved exactly (powers of two).
       Flooring a sub-ulp bound to one whole ulp of a coarse anchor
       destroys balls of the form q - epsilon whose error sits many
       limbs below their own mantissa ulp -- exactly what the
       cosine factors Q 2^E - A of the trigonometric bit-burst look
       like. */
    while (x->v < 1.0)
    {
        x->v = x->v * FBALL_D_B;
        x->a--;
    }
}

static fberr_t
fberr(double v, slong a)
{
    fberr_t e;
    e.v = v;
    e.a = a;
    fberr_reduce(&e);
    return e;
}

static fberr_t
fberr_add(fberr_t x, fberr_t y)
{
    fberr_t r;
    slong d;
    double t;

    if (x.v == 0.0) return y;
    if (y.v == 0.0) return x;

    if (x.a < y.a)
    {
        r = x; x = y; y = r;
    }

    d = x.a - y.a;
    t = (d > 40) ? 0.0 : ldexp(y.v, (int) (-FLINT_BITS * d));
    r.v = (x.v + t) * FBALL_EPS;
    r.a = x.a;
    fberr_reduce(&r);
    return r;
}

/* the error bound of x as an anchored pair */
static fberr_t
fberr_of(const fball_t x)
{
    return fberr(x->err, x->exp - x->size + x->erra);
}

/* strip high zero limbs (and low zero limbs when exact); zero
   mantissas collapse to size 0 keeping the anchor in exp */
static void
_fball_norm(fball_t x)
{
    while (x->size > 0 && x->d[x->size - 1] == 0)
    {
        x->size--;
        x->exp--;
    }

    if (x->size == 0)
    {
        x->negative = 0;
        return;
    }

    if (x->err == 0.0)
    {
        slong t;
        for (t = 0; x->d[t] == 0; t++)
            ;
        if (t > 0)
        {
            flint_mpn_copyi(x->d, x->d + t, x->size - t);
            x->size -= t;
        }
    }
}

/* Install the composed error bound e on x (whose mantissa, sign and
   exp/size are already set and top-normalized).  Low limbs are
   truncated until the bound, expressed in output ulps, is at most
   FBALL_ERR_MAX; each dropped run of limbs is worth < 1 ulp at the
   new anchor. */
static void
_fball_apply_err(fball_t x, fberr_t e)
{
    slong anc, drop, k;

    fberr_reduce(&e);

    if (e.v == 0.0)
    {
        x->err = 0.0;
        x->erra = 0;
        return;
    }

    anc = x->exp - x->size;
    k = e.a - anc;      /* e = v * B^k ulps of x */

    /* limbs to drop so that FLINT_BITS (k - drop) + log2(v)
       <= FBALL_ERR_MAX_EXP; frexp: v < 2^eb */
    {
        int eb;
        frexp(e.v, &eb);
        drop = k + (((slong) eb - FBALL_ERR_MAX_EXP)
            + FLINT_BITS - 1) / FLINT_BITS;
    }

    if (drop > 0)
    {
        double extra;

        if (drop > x->size - 1)
        {
            /* the bound swamps the mantissa: represent as
               0 +/- (mantissa magnitude + e); |x| < B^exp = one
               ulp at B^exp */
            e = fberr_add(e, fberr(1.0, x->exp));
            x->size = 0;
            x->exp = e.a;
            x->negative = 0;
            x->err = e.v;
            x->erra = 0;
            return;
        }

        flint_mpn_copyi(x->d, x->d + drop, x->size - drop);
        x->size -= drop;
        anc += drop;
        extra = 1.0;    /* dropped limbs: < 1 ulp at the new anchor */
        e = fberr_add(e, fberr(extra, anc));
        k = e.a - anc;
    }

    /* store the bound at its own anchor: no scaling, no underflow
       (a bound far below one ulp of the mantissa bottom -- a ball
       1 - epsilon -- keeps its full quality) */
    x->err = e.v;
    x->erra = e.a - anc;
    (void) k;

    if (x->size > 0 && x->d[x->size - 1] == 0)
        _fball_norm(x);
}

/* set x to the ball 0 +/- e */
static void
_fball_zero_err(fball_t x, fberr_t e)
{
    fberr_reduce(&e);
    x->size = 0;
    x->negative = 0;
    x->exp = e.a;
    x->err = e.v;
    x->erra = 0;
}

void
fball_init(fball_t x)
{
    x->d = NULL;
    x->alloc = 0;
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0.0;
    x->erra = 0;
}

void
fball_clear(fball_t x)
{
    flint_free(x->d);
}

void
fball_fit(fball_t x, slong k)
{
    if (x->alloc < k)
    {
        slong newalloc = FLINT_MAX(k, x->alloc + x->alloc / 2);
        x->d = flint_realloc(x->d, newalloc * sizeof(ulong));
        x->alloc = newalloc;
    }
}

void
fball_zero(fball_t x)
{
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0.0;
    x->erra = 0;
}

int
fball_is_zero_exact(const fball_t x)
{
    return x->size == 0 && x->err == 0.0;
}

void
fball_set_ui(fball_t x, ulong c)
{
    if (c == 0)
    {
        fball_zero(x);
        return;
    }

    fball_fit(x, 1);
    x->d[0] = c;
    x->size = 1;
    x->negative = 0;
    x->exp = 1;
    x->err = 0.0;
    x->erra = 0;
}

void
fball_set_si(fball_t x, slong c)
{
    fball_set_ui(x, (c >= 0) ? (ulong) c : -(ulong) c);
    x->negative = (c < 0);
}

void
fball_set(fball_t res, const fball_t x)
{
    if (res == x)
        return;
    fball_fit(res, x->size);
    flint_mpn_copyi(res->d, x->d, x->size);
    res->size = x->size;
    res->negative = x->negative;
    res->exp = x->exp;
    res->err = x->err;
    res->erra = x->erra;
}

void
fball_swap(fball_t x, fball_t y)
{
    FLINT_SWAP(fball_struct, *x, *y);
}

void
fball_add_error_ulps(fball_t x, double e)
{
    _fball_apply_err(x, fberr_add(fberr_of(x), fberr(e, x->exp - x->size)));
}

void
fball_add_error(fball_t x, double v, slong anc)
{
    _fball_apply_err(x, fberr_add(fberr_of(x), fberr(v, anc)));
}

void
fball_add_error_2exp_rel(fball_t x, slong e2)
{
    /* |value| < B^exp; add 2^e2 B^exp =
       2^(e2 mod FLINT_BITS) B^(exp + e2/FLINT_BITS), flooring the
       limb division so the bit remainder is >= 0 */
    slong q = e2 >> FBALL_LGB;
    int r = (int) (e2 - (q << FBALL_LGB));
    _fball_apply_err(x, fberr_add(fberr_of(x),
        fberr(ldexp(1.0, r), x->exp + q)));
}

/* ACCURATE LIMBS.  With errors normalized below 2^69 ulps, an
   inexact mantissa is accurate to all but its ~2 lowest limbs, so
   the mantissa length itself tracks accuracy; the target length of
   an output whose relative accuracy is limited by operand x is
   size + 1 when x is inexact (one guard limb; anything longer only
   computes noise) and unbounded when x is exact. */
/* accurate limbs of x: the mantissa length plus one, extended by
   however many whole limbs the error sits BELOW one ulp -- a ball
   like (q, err = 2^-800) arising from q - epsilon is accurate far
   beyond its short mantissa, and products against it must not be
   truncated to size + 1 limbs (rounding the extension down only
   ever underestimates the accuracy, which is safe) */
/* limb position of x's noise floor: the error anchor, extended
   downward by however many whole limbs a sub-ulp error bound sits
   below one ulp (mirroring _fball_acc; content above this line is
   meaningful and must not be discarded by additions) */
static slong
_fball_noise_floor(const fball_t x)
{
    int eb;
    frexp(x->err, &eb);     /* 1 <= err < 2^eb, eb in [1, 128] */
    return x->exp - x->size + x->erra + ((slong) eb - 1) / FLINT_BITS;
}

/* The pad keeps each operation's own truncation noise TWO limbs
   below the operand's noise floor instead of at it: a chain of m
   accuracy-capped operations then adds m 2^-(64 PAD) ulps of the
   floor to the radius rather than m whole ulps, and -- since
   interval radii over correlated errors (Karatsuba) grow faster
   than true errors -- the pad is also the number of limbs of TRUE
   accuracy that survive a radius overestimated by up to
   B^PAD -- without the pad, a
   ball imported with a 96-ulp radius (a series-leaf trig factor)
   drifts to several hundred ulps across one bit-burst
   accumulation. */
#define FBALL_ACC_PAD 6

static slong
_fball_acc(const fball_t x, slong n)
{
    if (x->err == 0.0)
        return WORD_MAX;
    return x->exp - _fball_noise_floor(x) + 1 + FBALL_ACC_PAD;
}

/* res = a * b truncated to about n limbs.  Full products are used
   when at most a few limbs longer than the kept part (the low limbs
   then round down into <= 1 ulp), otherwise flint_mpn_mulmid computes
   just the kept window plus one guard limb below, whose one-sided
   deficit of up to min(an, bn) + 2 ulps one limb above the window
   bottom is absorbed into the radius (and immediately truncated away
   by the error normalization). */
void
fball_mul(fball_t res, const fball_t a, const fball_t b, slong n)
{
    slong sa = a->size, sb = b->size, ps, keep;
    slong rexp;
    int rneg;
    fberr_t e;
    TMP_INIT;

    /* cross terms |a| errb + |b| erra + erra errb, valid also for
       zero mantissas.  The magnitudes are bounded through the top
       limb, |a| < (a_top + 1) B^(exp-1), rather than by B^exp:
       the blanket bound overstates by up to a full limb whenever
       the top limb is small (maximally for a multiplier equal to
       1), and that factor-B pessimism compounds across a chain of
       multiplications.  The double conversion of the top limb can
       round down by up to 2^11, covered by the (1 + 2^-50) fudge
       before the + 1. */
    e = fberr(0.0, 0);
    if (b->err != 0.0)
    {
        double ma = (sa == 0) ? 0.0
            : (double) a->d[sa - 1] * (1.0 + 0x1p-50) + 1.0;
        if (sa == 0)
            e = fberr_add(e, fberr(b->err,
                    a->exp + b->exp - b->size + b->erra));
        else
            e = fberr_add(e, fberr(b->err * ma,
                    (a->exp - 1) + b->exp - b->size + b->erra));
    }
    if (a->err != 0.0)
    {
        double mb = (sb == 0) ? 0.0
            : (double) b->d[sb - 1] * (1.0 + 0x1p-50) + 1.0;
        if (sb == 0)
            e = fberr_add(e, fberr(a->err,
                    b->exp + a->exp - a->size + a->erra));
        else
            e = fberr_add(e, fberr(a->err * mb,
                    (b->exp - 1) + a->exp - a->size + a->erra));
        if (b->err != 0.0)
            e = fberr_add(e, fberr(a->err * b->err,
                    a->exp - a->size + a->erra
                    + b->exp - b->size + b->erra));
    }

    if (sa == 0 || sb == 0)
    {
        _fball_zero_err(res, e);
        return;
    }

    ps = sa + sb;
    keep = FLINT_MIN(n, ps);
    keep = FLINT_MIN(keep, _fball_acc(a, n));
    keep = FLINT_MIN(keep, _fball_acc(b, n));
    keep = FLINT_MAX(keep, 1);

    rexp = a->exp + b->exp;
    rneg = a->negative ^ b->negative;

    if (ps <= keep + 4 || ps <= keep + keep / 16)
    {
        /* full product, written straight into the destination
           buffer unless it would overwrite an operand */
        slong drop = ps - keep;
        int aliased = (res == a || res == b);
        nn_ptr t;
        TMP_START;
        if (aliased)
            t = TMP_ALLOC(ps * sizeof(ulong));
        else
        {
            fball_fit(res, ps);
            t = res->d;
        }

        if (a->d == b->d && sa == sb)
            flint_mpn_sqr(t, a->d, sa);
        else if (sa >= sb)
            flint_mpn_mul(t, a->d, sa, b->d, sb);
        else
            flint_mpn_mul(t, b->d, sb, a->d, sa);

        if (drop > 0)
        {
            /* keep exactness when the dropped tail is zero */
            slong i;
            int inexact = 0;
            for (i = 0; i < drop && !inexact; i++)
                inexact = (t[i] != 0);
            if (inexact)
                e = fberr_add(e, fberr(1.0, rexp - ps + drop));
        }

        if (aliased)
        {
            fball_fit(res, keep);
            flint_mpn_copyi(res->d, t + drop, keep);
        }
        else if (drop > 0)
            flint_mpn_copyi(res->d, res->d + drop, keep);
        res->size = keep;
        res->exp = rexp;
        res->negative = rneg;
        res->err = e.v;     /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_err(res, e);
        TMP_END;
    }
    else
    {
        /* windowed middle product: window [ps - keep - 1, ps) */
        slong zlo = ps - keep - 1, zn = keep + 1;
        int aliased = (res == a || res == b);
        nn_ptr t;
        TMP_START;
        if (aliased)
            t = TMP_ALLOC(zn * sizeof(ulong));
        else
        {
            fball_fit(res, zn);
            t = res->d;
        }

        flint_mpn_mulmid(t, a->d, sa, b->d, sb, zlo, ps);

        /* deficit (below the true window) plus the discarded low
           limbs: < (min + 2) B + 1 ulps of the window bottom */
        e = fberr_add(e, fberr(
            (double) (FLINT_MIN(FLINT_MIN(sa, sb), zlo) + 2),
            rexp - ps + zlo + 1));
        e = fberr_add(e, fberr(1.0, rexp - ps + zlo));

        if (aliased)
        {
            fball_fit(res, zn);
            flint_mpn_copyi(res->d, t, zn);
        }
        res->size = zn;
        res->exp = rexp;
        res->negative = rneg;
        res->err = e.v;     /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_err(res, e);
        TMP_END;
    }
}

void
fball_mul_ui(fball_t res, const fball_t a, ulong c, slong n)
{
    fberr_t e;
    nn_ptr t;
    ulong cy;
    slong sa = a->size, keep, drop;

    if (c == 0)
    {
        fball_zero(res);
        return;
    }

    e = fberr(a->err * (double) c * FBALL_EPS,
        a->exp - a->size + a->erra);

    if (sa == 0)
    {
        _fball_zero_err(res, fberr_add(e, fberr(0.0, 0)));
        return;
    }

    /* works in the destination buffer, in place when aliased
       (fball_fit preserves contents across reallocation) */
    fball_fit(res, sa + 1);
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
                e = fberr_add(e, fberr(1.0, anc + drop));
            flint_mpn_copyi(t, t + drop, keep);
        }

        res->size = keep;
        res->exp = anc + ts;
        res->negative = a->negative;
        res->err = e.v;     /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_err(res, e);
    }
}

/* signed addition core: res = (-1)^aneg |a| + (-1)^bneg |b| */
static void
_fball_add(fball_t res, const fball_t a, int aneg, const fball_t b,
    int bneg, slong n)
{
    slong top, bot, span, oa, ob;
    fberr_t e;
    nn_ptr t;
    int rneg, trunc_a = 0, trunc_b = 0;
    TMP_INIT;

    e = fberr_add(fberr_of(a), fberr_of(b));

    if (a->size == 0 && b->size == 0)
    {
        _fball_zero_err(res, e);
        return;
    }

    if (a->size == 0 || b->size == 0)
    {
        /* e already contains both radii; apply_err overwrites err */
        if (a->size == 0)
        {
            fball_set(res, b);
            res->negative = bneg;
        }
        else
        {
            fball_set(res, a);
            res->negative = aneg;
        }
        _fball_apply_err(res, e);
        return;
    }

    top = FLINT_MAX(a->exp, b->exp) + 1;
    bot = FLINT_MIN(a->exp - a->size, b->exp - b->size);
    bot = FLINT_MAX(bot, top - (n + 2));
    /* limbs well below either operand's noise floor are noise
       (keep FBALL_ACC_PAD limbs under the floor so this truncation
       stays negligible against the operand's radius) */
    if (a->err != 0.0)
        bot = FLINT_MAX(bot, _fball_noise_floor(a) - FBALL_ACC_PAD);
    if (b->err != 0.0)
        bot = FLINT_MAX(bot, _fball_noise_floor(b) - FBALL_ACC_PAD);

    span = top - bot;

    /* The |a| window is built directly in the destination buffer:
       a temporary is needed only when res aliases b (whose limbs
       would be overwritten before the combine); when res aliases a
       the mantissa is shifted within its own buffer.  Only the
       gaps around |a| are zeroed, not the whole span. */
    oa = (a->exp - a->size) - bot;

    TMP_START;
    if (res == b)
        t = TMP_ALLOC(span * sizeof(ulong));
    else
    {
        fball_fit(res, span);
        t = res->d;
    }

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
            flint_mpn_zero(t, oa);
        flint_mpn_zero(t + oa + a->size, span - oa - a->size);
    }
    else if (a->size + oa > 0)
    {
        /* low limbs of a truncated; copying downward is safe in
           place (ascending copy, dst below src) */
        flint_mpn_copyi(t, (res == a ? t : a->d) - oa, a->size + oa);
        flint_mpn_zero(t + a->size + oa, span - a->size - oa);
        trunc_a = 1;
    }
    else
    {
        trunc_a = (a->size > 0);
        flint_mpn_zero(t, span);
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
                ulong cy = mpn_add(t + ob, t + ob, span - ob, bd, bs);
                FLINT_ASSERT(cy == 0);  /* top slot reserved */
                (void) cy;
            }
            else
            {
                ulong bw = mpn_sub(t + ob, t + ob, span - ob, bd, bs);
                if (bw)
                {
                    /* |b| window > |a| window: negate */
                    mpn_neg(t, t, span);
                    rneg = bneg;
                }
            }
        }
        else
        {
            trunc_b |= (b->size > 0);
        }
    }

    if (trunc_a)
        e = fberr_add(e, fberr(1.0, bot));
    if (trunc_b)
        e = fberr_add(e, fberr(1.0, bot));

    if (res == b)
    {
        fball_fit(res, span);
        flint_mpn_copyi(res->d, t, span);
    }
    res->size = span;
    res->exp = top;
    res->negative = rneg;
    res->err = e.v;     /* norm must know inexactness */
    _fball_norm(res);
    if (res->size == 0)
        _fball_zero_err(res, fberr_add(e, fberr(0.0, 0)));
    else
        _fball_apply_err(res, e);
    TMP_END;
}

void
fball_add(fball_t res, const fball_t a, const fball_t b, slong n)
{
    _fball_add(res, a, a->negative, b, b->negative, n);
}

void
fball_sub(fball_t res, const fball_t a, const fball_t b, slong n)
{
    _fball_add(res, a, a->negative, b, !b->negative, n);
}

/* rigorous OVERestimate of the relative radius err / |x| of a ball
   with nonzero mantissa: |x| >= B^(exp - 1), so
   rel <= err B^(erra + 1 - size); the ldexp window is clamped on
   the small side to a value exceeding anything the invariant
   err < 2^128 allows there, keeping the estimate one-sided */
static double
_fball_rel_bound(const fball_t x)
{
    slong k = x->erra + 1 - x->size;

    if (x->err == 0.0)
        return 0.0;
    if (FLINT_BITS * k > 900)
        return HUGE_VAL;
    if (FLINT_BITS * k < -900)
        return 0x1p-700;    /* true bound < 2^(128 - 900) */
    return ldexp(x->err, (int) (FLINT_BITS * k)) * FBALL_EPS;
}

/* res = a / b.  Both mantissas are read as fractions and shifted
   left, together, by the divisor's leading-zero count -- the same
   quotient with a top-bit-normalized divisor, so the Newton error
   4 B^-n / d is at most 8 ulps (the numerator gains a zero top limb
   to stay below 1).  Operand errors enter as
   err_a / |b| + |a| err_b / |b|^2, with |b| bounded below by
   B^(b.exp - 1) and the -err_b correction of the denominator covered
   by a 1.01 factor, which is rigorous because the relative radius of
   b is CHECKED to be below 2^-30 (a divisor ball this wide is a
   usage error: the mantissa would be pure noise). */
void
fball_div(fball_t res, const fball_t a, const fball_t b, slong n)
{
    slong sa = a->size, sb = b->size, nd, qn;
    slong rexp;
    int rneg, c;
    fberr_t e;
    nn_ptr num, den, q;
    TMP_INIT;

    if (sb == 0)
        flint_throw(FLINT_ERROR, "fball_div: division by zero ball\n");
    if (_fball_rel_bound(b) > 0x1p-30)
        flint_throw(FLINT_ERROR,
            "fball_div: divisor relative radius above 2^-30\n");

    e = fberr(0.0, 0);
    if (a->err != 0.0)
        e = fberr_add(e, fberr(a->err * 1.01,
                (a->exp - sa) + a->erra + 1 - b->exp));
    if (b->err != 0.0)
        e = fberr_add(e, fberr(b->err * 1.01,
                a->exp + (b->exp - sb) + b->erra + 2 - 2 * b->exp));

    if (sa == 0)
    {
        _fball_zero_err(res, e);
        return;
    }

    nd = FLINT_MIN(n, FLINT_MIN(_fball_acc(a, n), _fball_acc(b, n)));
    nd = FLINT_MAX(nd, 2);
    qn = nd + 2;

    rneg = a->negative ^ b->negative;

    /* the quotient goes straight into the destination buffer when
       it does not alias an operand; shift copies of the operands
       are made only when the divisor needs top-bit normalization */
    c = flint_clz(b->d[sb - 1]);

    TMP_START;
    {
        int aliased = (res == a || res == b);
        slong tn = (c ? (sa + 1 + sb) : 0) + (aliased ? qn : 0);
        nn_ptr t = (tn > 0) ? TMP_ALLOC(tn * sizeof(ulong)) : NULL;

        if (!aliased)
        {
            fball_fit(res, qn);
            q = res->d;
        }
        else
            q = t + (c ? sa + 1 + sb : 0);

        if (c == 0)
        {
            /* divisor already top-bit normalized:
               quotient = (ma B^-sa) / (mb B^-sb) */
            fixed_div_newton(q, a->d, sa, b->d, sb, nd);
            rexp = a->exp - b->exp + 2;
        }
        else
        {
            /* shift numerator and denominator together (the
               numerator gains a zero top limb to remain a fraction
               below 1); quotient = (ma 2^c B^-(sa+1)) /
               (mb 2^c B^-sb), one limb smaller than above */
            num = t;
            den = t + sa + 1;
            num[sa] = mpn_lshift(num, a->d, sa, c);
            mpn_lshift(den, b->d, sb, c);
            fixed_div_newton(q, num, sa + 1, den, sb, nd);
            rexp = a->exp - b->exp + 3;
        }

        /* Newton error <= 4 B^-nd / den <= 8 ulps of
           B^(rexp - qn) */
        e = fberr_add(e, fberr(8.0, rexp - qn));

        if (aliased)
        {
            fball_fit(res, qn);
            flint_mpn_copyi(res->d, q, qn);
        }
    }
    res->size = qn;
    res->exp = rexp;
    res->negative = rneg;
    res->err = e.v;     /* norm must know inexactness */
    _fball_norm(res);
    if (res->size == 0)
        _fball_zero_err(res, fberr_add(e, fberr(0.0, 0)));
    else
        _fball_apply_err(res, e);
    TMP_END;
}

void
fball_rsqrt_ui(fball_t res, ulong c, slong n)
{
    FLINT_ASSERT(c >= 2);

    fball_fit(res, n);
    fixed_rsqrt_ui_newton(res->d, c, n);
    res->size = n;
    res->exp = 0;
    res->negative = 0;
    res->err = 1.0;     /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_err(res, fberr(2.0, -n));
}

void
fball_set_mpn_2exp(fball_t x, nn_srcptr p, slong len, slong ebits)
{
    slong b, q;

    while (len > 0 && p[len - 1] == 0)
        len--;
    if (len == 0)
    {
        /* keep the frame: a zero import anchors its ulp at
           B^ceil(ebits / FLINT_BITS) >= 2^ebits, so that a radius
           subsequently attached in ulps (a truncated-to-zero
           quantity known to |value| < k 2^ebits) lands at the
           intended scale instead of at B^0 */
        fball_zero(x);
        x->exp = ebits / FLINT_BITS + ((ebits % FLINT_BITS) > 0);
        return;
    }

    b = ebits % FLINT_BITS;
    if (b < 0)
        b += FLINT_BITS;
    q = (ebits - b) / FLINT_BITS;

    fball_fit(x, len + 1);
    if (b)
    {
        x->d[len] = mpn_lshift(x->d, p, len, (int) b);
        x->size = len + (x->d[len] != 0);
    }
    else
    {
        flint_mpn_copyi(x->d, p, len);
        x->size = len;
    }
    x->negative = 0;
    x->exp = x->size + q;
    x->err = 0.0;
    x->erra = 0;
    _fball_norm(x);
}

double
fball_get_fixed(nn_ptr y, slong wn, const fball_t x)
{
    double e;
    slong sh;

    /* radius in output ulps: err B^(exp - size + wn), saturating.
       Both saturation directions must OVERestimate: far above, the
       bound is useless anyway (HUGE_VAL); far below, err < 2^128
       at k <= -k0 limbs is at most 2^(128 - FLINT_BITS k0), so the
       clamp constant must exceed that, while the ldexp in the live
       window (|FLINT_BITS k| <= 900) cannot underflow */
    if (x->err == 0.0)
        e = 0.0;
    else
    {
        slong k = x->exp - x->size + x->erra + wn;
        if (FLINT_BITS * k > 900)
            e = HUGE_VAL;
        else if (FLINT_BITS * k < -900)
            e = 0x1p-700;   /* true bound < 2^(128 - 900) = 2^-772 */
        else
            e = ldexp(x->err, (int) (FLINT_BITS * k)) * FBALL_EPS;
    }

    if (x->size == 0)
    {
        flint_mpn_zero(y, wn);
        return e + 1.0;
    }

    if (x->negative)
    {
        /* nonpositive point value in [0,1) ball: clamp; the true
           value is within (radius - value) <= radius of 0 */
        flint_mpn_zero(y, wn);
        return e;
    }

    FLINT_ASSERT(x->exp <= 0);  /* value < 1 */

    sh = wn + x->exp - x->size;
    flint_mpn_zero(y, wn);
    if (sh >= 0)
    {
        FLINT_ASSERT(sh + x->size <= wn);
        flint_mpn_copyi(y + sh, x->d, x->size);
        return e;
    }
    /* drop the low -sh limbs: one-sided < 1 output ulp */
    if (-sh < x->size)
        flint_mpn_copyi(y, x->d - sh, x->size + sh);
    return e + 1.0;
}

/* If the ball x, known to satisfy 0 <= x < 1, determines
   floor(x B^n) uniquely, write it to (y, n) and return 1; otherwise
   return 0 (caller retries at higher precision).  Everything is
   limb-aligned: the mantissa splits at the output grid into the kept
   top part M and a t-limb tail L, and the radius -- overestimated by
   the pure power of two 2^s, s = eb + FLINT_BITS erra with
   err < 2^eb from frexp -- must fit strictly inside the tail on both
   sides: 2^s <= L (the true value cannot borrow below the grid line)
   and L < B^t - 2^s (nor carry above the next), both checked by bit
   scans without materializing L + 2^s. */
int
fball_get_fixed_floor(nn_ptr y, slong n, const fball_t x)
{
    slong sh, t, s, i, hi;
    int eb;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(x->exp <= 0);
    FLINT_ASSERT(x->size == 0 || !x->negative);

    flint_mpn_zero(y, n);

    if (x->size == 0)
    {
        /* 0 +/- err: the floor is 0 iff the radius stays strictly
           below B^-n: err B^(exp + erra) < B^-n */
        if (x->err == 0.0)
            return 1;
        frexp(x->err, &eb);
        return (slong) eb + FLINT_BITS * (x->exp + x->erra + n) <= 0;
    }

    sh = x->exp - x->size + n;  /* x B^n = (d, size) B^sh */
    t = -sh;

    if (t <= 0)
    {
        /* mantissa entirely on or above the grid: exact values
           only (any radius straddles the grid line the value sits
           on) */
        if (x->err != 0.0)
            return 0;
        FLINT_ASSERT(x->size + sh <= n);
        flint_mpn_copyi(y + sh, x->d, x->size);
        return 1;
    }

    if (t < x->size)
        flint_mpn_copyi(y, x->d + t, x->size - t);
    /* else the whole value sits below the grid: floor 0, and the
       t-limb tail is the mantissa padded above with zero limbs */

    if (x->err == 0.0)
        return 1;

    frexp(x->err, &eb);
    s = (slong) eb + FLINT_BITS * x->erra;
    s = FLINT_MAX(s, 0);

    if (s >= FLINT_BITS * t - 1)
        return 0;

    /* the tail L, read as t limbs with zeros above the mantissa */
#define FBALL_TAIL_LIMB(I_) ((I_) < x->size ? x->d[I_] : UWORD(0))

    /* low side: L >= 2^s, i.e. some bit at position >= s is set */
    hi = -1;
    for (i = t - 1; i >= 0; i--)
        if (FBALL_TAIL_LIMB(i) != 0)
        {
            hi = i;
            break;
        }
    if (hi < 0 || FLINT_BITS * hi
            + (slong) FLINT_BIT_COUNT(FBALL_TAIL_LIMB(hi)) <= s)
        return 0;

    /* high side: L < B^t - 2^s, i.e. some bit at position >= s is
       clear */
    {
        slong ls = s / FLINT_BITS;
        int sb = (int) (s % FLINT_BITS);

        for (i = t - 1; i > ls; i--)
            if (FBALL_TAIL_LIMB(i) != UWORD_MAX)
                return 1;
        return (FBALL_TAIL_LIMB(ls) >> sb) != (UWORD_MAX >> sb);
    }
#undef FBALL_TAIL_LIMB
}

void
fball_get_arb(arb_t res, const fball_t x)
{
    fmpz_t m;
    fmpz_init(m);

    if (x->size == 0)
        arf_zero(arb_midref(res));
    else
    {
        fmpz_set_ui_array(m, x->d, x->size);
        if (x->negative)
            fmpz_neg(m, m);
        arf_set_fmpz(arb_midref(res), m);
        arf_mul_2exp_si(arb_midref(res), arb_midref(res),
            FLINT_BITS * (x->exp - x->size));
    }

    if (x->err == 0.0)
        mag_zero(arb_radref(res));
    else
    {
        mag_set_d(arb_radref(res), x->err);
        mag_mul_2exp_si(arb_radref(res), arb_radref(res),
            FLINT_BITS * (x->exp - x->size + x->erra));
    }

    fmpz_clear(m);
}

void
fball_print(const fball_t x)
{
    arb_t t;
    arb_init(t);
    fball_get_arb(t, x);
    flint_printf("[size %wd, exp %wd, err %g @%wd] = ",
        x->size, x->exp, x->err, x->erra);
    arb_printd(t, 20);
    flint_printf("\n");
    arb_clear(t);
}

/* value *= 2^e (bit-level shift; everything else in fball is
   limb-aligned) */
void
fball_mul_2exp_si(fball_t x, slong e)
{
    slong q = e >> FBALL_LGB;             /* floor limb division */
    int r = (int) (e - (q << FBALL_LGB)); /* 0 <= r < FLINT_BITS */

    if (x->size == 0)
    {
        x->exp += q + (r != 0);
        /* the radius scales exactly; anchor moved by q(+1) limbs */
        if (x->err != 0.0)
            _fball_apply_err(x,
                fberr(ldexp(x->err, r - ((r != 0) << FBALL_LGB)),
                x->exp + x->erra));
        return;
    }

    x->exp += q;
    if (r != 0)
    {
        ulong cy;
        double err = x->err;
        slong erra0 = x->erra;

        fball_fit(x, x->size + 1);
        cy = mpn_lshift(x->d, x->d, x->size, r);
        if (cy)
        {
            x->d[x->size] = cy;
            x->size++;
            x->exp++;
        }
        _fball_norm(x);
        if (err != 0.0)
            _fball_apply_err(x, fberr(ldexp(err, r),
                x->exp - x->size + erra0));
        else
        {
            x->err = 0.0;
            x->erra = 0;
        }
    }
}

/* shared alignment for sqrt/rsqrt: express x = ahat * B^E with E
   even, ahat in [B^-2, 1) given as an an-limb fraction (aa, an);
   the input radius keeps its anchor: err ulps of B^(E - an) */
#define FBALL_SQRT_ALIGN                                        \
    slong E = x->exp, an = x->size, nd;                         \
    nn_ptr aa;                                                  \
    fberr_t e;                                                  \
    TMP_INIT;                                                   \
    FLINT_ASSERT((const void *) res != (const void *) x);       \
    FLINT_ASSERT(x->size > 0 && !x->negative);                  \
    if (_fball_rel_bound(x) > 0x1p-30)                          \
        flint_throw(FLINT_ERROR,                                \
            "fball sqrt: relative radius above 2^-30\n");       \
    nd = FLINT_MIN(n, _fball_acc(x, n));                        \
    nd = FLINT_MAX(nd, 2);                                      \
    TMP_START;                                                  \
    if (E & 1)                                                  \
    {                                                           \
        aa = TMP_ALLOC((an + 1) * sizeof(ulong));               \
        flint_mpn_copyi(aa, x->d, an);                          \
        aa[an] = 0;                                             \
        an++;                                                   \
        E++;                                                    \
    }                                                           \
    else                                                        \
        aa = (nn_ptr) x->d;

/* res = 1/sqrt(x); x > 0, n-limb target */
void
fball_rsqrt(fball_t res, const fball_t x, slong n)
{
    FBALL_SQRT_ALIGN

    /* operand error through the derivative: with
       x >= B^(E-1-q), q = 1 iff the alignment padded an odd
       exponent (zero top limb), |d(1/sqrt x)/dx| = 1/(2 x^(3/2))
       <= B^(3(1+q-E)/2)/2, so
       |Delta| <= 1.02 err sqrt(B) B^(1+q-an-E/2+erra) B^(q/2).
       Anchoring at the blanket magnitude bound instead loses ~1.5
       limbs of claimed accuracy PER CALL, which the noise-floor
       truncation of subsequent operations turns into real mantissa
       loss (fatal for AGM chains of dozens of square roots). */
    e = fberr(0.0, 0);
    if (x->err != 0.0)
    {
        int q = (int) (x->exp & 1);
        e = fberr(x->err * 1.02 * (q ? FBALL_D_B : FBALL_D_SQRTB),
            1 + q - an - E / 2 + x->erra);
    }

    fball_fit(res, nd + 2);
    fixed_rsqrt_newton(res->d, aa, an, nd);

    /* Newton: <= 4 B^-nd / sqrt(ahat) <= 4 B^(1 - nd), scaled
       B^(-E/2) */
    e = fberr_add(e, fberr(4.0, 1 - nd - E / 2));

    res->size = nd + 2;
    res->exp = 2 - E / 2;
    res->negative = 0;
    res->err = e.v;     /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_err(res, e);
    TMP_END;
}

/* res = sqrt(x); x > 0, n-limb target */
void
fball_sqrt(fball_t res, const fball_t x, slong n)
{
    FBALL_SQRT_ALIGN

    /* operand error through the derivative: with x >= B^(E-1-q)
       as above, |d(sqrt x)/dx| = 1/(2 sqrt x) <= B^((1+q-E)/2)/2:
       |Delta| <= 1.02 err B^((q - 1)/2) B^(1-an+E/2+erra) */
    e = fberr(0.0, 0);
    if (x->err != 0.0)
    {
        int q = (int) (x->exp & 1);
        e = fberr(x->err * 1.02 * (q ? 1.0 : FBALL_D_SQRTBINV),
            1 - an + E / 2 + x->erra);
    }

    fball_fit(res, nd + 2);
    fixed_sqrt_newton(res->d, aa, an, nd);

    /* Newton: <= 4 B^-nd / sqrt(ahat) <= 4 B^(1 - nd), scaled
       B^(E/2) */
    e = fberr_add(e, fberr(4.0, 1 - nd + E / 2));

    res->size = nd + 2;
    res->exp = 2 + E / 2;
    res->negative = 0;
    res->err = e.v;     /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_err(res, e);
    TMP_END;
}

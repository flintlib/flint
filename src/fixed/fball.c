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
#include "fixed.h"
#include "arb.h"

/* ERROR BOUNDS.  A ball's radius is an integer count err < B of ulps
   of its bottom limb.  During the composition of a bound the
   contributions (the radii of the operands scaled by magnitudes, the
   truncations) live at limb anchors that may differ, so they are
   carried as pairs

       e = (hi B + lo) B^a,   a 128-bit count at the limb anchor a,

   and combined rounding UP in limb arithmetic: a sum aligns the lower
   anchor to the higher one, dividing its count by B^d and rounding up
   (one unit when d >= 2), so that every term costs at most one unit
   of the higher anchor beyond its value; a product of a radius by a
   magnitude (a top limb plus one) is one umul_ppmm.  Installing a
   bound on a ball then converts it to units of the mantissa's bottom
   limb, truncating limbs when it reaches B units (the dropped material
   is worth one more unit) and padding zero limbs when it lies below
   one ulp.  The bounds of the divisions and roots are composed in
   doubles (fberr_t, anchored the same way, rounding up with the fudge
   factor FBALL_EPS) and converted at the end. */

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

/* a 128-bit count at a limb anchor */
typedef struct { ulong hi, lo; slong a; } fbnd_t;

static const fbnd_t fbnd_zero = { 0, 0, 0 };

FLINT_FORCE_INLINE int
fbnd_is_zero(fbnd_t x)
{
    return (x.hi | x.lo) == 0;
}

/* the radius of x */
FLINT_FORCE_INLINE fbnd_t
fbnd_of(const fball_t x)
{
    fbnd_t e;
    e.hi = 0;
    e.lo = x->err;
    e.a = (x->size == 0) ? x->exp : x->exp - x->size;
    return e;
}

/* (hi, lo) B^a */
FLINT_FORCE_INLINE fbnd_t
fbnd(ulong hi, ulong lo, slong a)
{
    fbnd_t e;
    e.hi = hi;
    e.lo = lo;
    e.a = a;
    return e;
}

/* MAGNITUDES.  The magnitude of a ball with nonzero mantissa is
   bounded through the leading FLINT_BITS bits m of its mantissa (across
   a limb boundary), normalized with the top bit set by a shift s:

       m 2^-s B^(exp - 1) <= |x| < (m + 1) 2^-s B^(exp - 1),

   so the bounds overstate by at most 2^-63 relatively.  Bounding
   through the top limb alone, |x| < (top + 1) B^(exp - 1), overstates
   by a factor up to 2 for a top limb of 1 (every value in [1, 2),
   the working set of any iteration near 1) and by up to B in the
   other direction; with the radius held as a whole count of ulps
   that pessimism turns directly into truncated limbs. */
FLINT_FORCE_INLINE ulong
_fball_lead(const fball_t x, int * s)
{
    ulong m = x->d[x->size - 1];
    int c = flint_clz(m);
    if (c != 0)
    {
        m <<= c;
        if (x->size > 1)
            m |= x->d[x->size - 2] >> (FLINT_BITS - c);
    }
    *s = c;
    return m;
}

/* lower and upper bounds of |x| / B^(exp - 1) as doubles in [1, B] */
FLINT_FORCE_INLINE double
_fball_mag_lo(const fball_t x)
{
    int s;
    ulong m = _fball_lead(x, &s);
    return ldexp((double) m * (1.0 - 0x1p-52), -s);
}

FLINT_FORCE_INLINE double
_fball_mag_hi(const fball_t x)
{
    int s;
    ulong m = _fball_lead(x, &s);
    return ldexp((double) m * (1.0 + 0x1p-52) + 1.0, -s);
}

/* v |x| B^(a + 1 - exp) bounded: v (m + 1) 2^-s B^a, rounded up */
FLINT_FORCE_INLINE fbnd_t
fbnd_mul_mag(ulong v, const fball_t x, slong a)
{
    fbnd_t e;
    int s;
    ulong m = _fball_lead(x, &s), r;

    umul_ppmm(e.hi, e.lo, v, m);
    add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), v);
    if (s != 0)
    {
        r = e.lo & ((UWORD(1) << s) - 1);
        e.lo = (e.lo >> s) | (e.hi << (FLINT_BITS - s));
        e.hi >>= s;
        if (r != 0)
            add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(1));
    }
    e.a = a;
    return e;
}

/* x + y rounded up: the lower anchor is aligned to the higher one, its
   count divided by B^d rounding up (one unit for d >= 2); a carry out
   of 128 bits moves the anchor up a limb */
FLINT_FORCE_INLINE fbnd_t
fbnd_add(fbnd_t x, fbnd_t y)
{
    fbnd_t r;
    slong d;
    ulong yhi, ylo, cy;

    if (fbnd_is_zero(x)) return y;
    if (fbnd_is_zero(y)) return x;

    if (x.a < y.a)
    {
        r = x; x = y; y = r;
    }
    d = x.a - y.a;
    if (d == 0)
    {
        yhi = y.hi;
        ylo = y.lo;
    }
    else if (d == 1)
    {
        ylo = y.hi + (y.lo != 0);
        yhi = (ylo == 0);            /* y.hi = B - 1 with a nonzero low */
    }
    else
    {
        ylo = 1;
        yhi = 0;
    }
    add_ssaaaa(r.hi, r.lo, x.hi, x.lo, yhi, ylo);
    cy = (r.hi < x.hi) || (r.hi == x.hi && r.lo < x.lo);
    r.a = x.a;
    if (cy)
    {
        /* the sum is 2^128 + (r.hi, r.lo): rounded up by a limb */
        ylo = r.hi + (r.lo != 0);
        r.hi = 1 + (ylo < r.hi);
        r.lo = ylo;
        r.a++;
    }
    return r;
}

/* the count of e as a single limb, rounded up, with its anchor: e is
   below 2^129 by construction, so at most two limbs move */
FLINT_FORCE_INLINE ulong
fbnd_reduce(fbnd_t * e)
{
    if (e->hi != 0)
    {
        ulong v = e->hi + (e->lo != 0);
        e->a++;
        if (v == 0)
        {
            v = 1;
            e->a++;
            /* if e->hi + 1 overflowed, e was in [2^128 - B + 1, 2^128):
               one unit of B^(a + 2) covers it */
        }
        e->hi = 0;
        e->lo = v;
    }
    return e->lo;
}

/* --- the double-anchored bounds of the divisions and roots --- */

typedef struct { double v; slong a; } fberr_t;

/* v == 0 or 1 <= v < 2^128, restored by rebasing the anchor (exact
   power-of-two scalings; a limb pushed below the retained part is
   rounded up to one unit) */
static void
fberr_reduce(fberr_t * x)
{
    if (x->v >= 0x1p128)
    {
        x->v = x->v * FBALL_D_BINV + 1.0;
        x->a++;
        while (x->v >= 0x1p128)
        {
            x->v = x->v * FBALL_D_BINV + 1.0;
            x->a++;
        }
    }
    else if (x->v < 1.0 && x->v != 0.0)
    {
        do
        {
            x->v = x->v * FBALL_D_B;
            x->a--;
        } while (x->v < 1.0);
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

/* x + y, rounded up: aligned to the higher anchor, the lower operand
   scaled by B^-d for d = 1, 2 (exact) and dropped for d >= 3, where
   it is below 2^-64 of the higher one and thus inside the fudge
   factor FBALL_EPS applied to the sum */
static fberr_t
fberr_add(fberr_t x, fberr_t y)
{
    fberr_t r;
    slong d;

    if (x.v == 0.0) return y;
    if (y.v == 0.0) return x;

    if (x.a < y.a)
    {
        r = x; x = y; y = r;
    }

    d = x.a - y.a;
    if (d == 0)
        r.v = x.v + y.v;
    else if (d == 1)
        r.v = x.v + y.v * FBALL_D_BINV;
    else if (d == 2)
        r.v = x.v + y.v * (FBALL_D_BINV * FBALL_D_BINV);
    else
        r.v = x.v;
    r.v *= FBALL_EPS;
    r.a = x.a;
    if (r.v >= 0x1p128)
    {
        r.v = r.v * FBALL_D_BINV + 1.0;
        r.a++;
    }
    return r;
}

/* a double-anchored bound as a limb count: v < 2^128 is split at
   B, the low part rounded up */
static fbnd_t
fbnd_of_fberr(fberr_t e)
{
    fbnd_t r;
    double hi, lo;

    fberr_reduce(&e);
    if (e.v == 0.0)
        return fbnd_zero;
    hi = floor(e.v * FBALL_D_BINV);
    lo = ceil(e.v - hi * FBALL_D_B);
    r.a = e.a;
    r.hi = (ulong) hi;
    if (lo >= FBALL_D_B)
    {
        r.lo = 0;
        r.hi++;
        if (r.hi == 0)
        {
            /* 2^128: one unit two limbs up */
            r.lo = 1;
            r.a += 2;
        }
    }
    else
        r.lo = (ulong) lo;
    return r;
}

/* strip high zero limbs (and low zero limbs when exact); zero
   mantissas collapse to size 0 keeping the anchor in exp */
static void
_fball_norm_slow(fball_t x)
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

    if (x->err == 0)
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

FLINT_FORCE_INLINE void
_fball_norm(fball_t x)
{
    /* the common case: nonzero top limb, and a nonzero bottom limb
       unless inexact (where low zero limbs are kept as padding) */
    if (x->size > 0 && x->d[x->size - 1] != 0 && (x->err != 0 || x->d[0] != 0))
        return;
    _fball_norm_slow(x);
}

/* Install the composed bound e on x (whose mantissa, sign and
   exp/size are already set and top-normalized), as a count of ulps of
   the bottom limb: a count of B or more truncates the mantissa by as
   many limbs as needed (each dropped run of limbs is worth one unit),
   a bound below one ulp pads the mantissa with zero limbs down to its
   scale. */
static void
_fball_apply_bnd_slow(fball_t x, fbnd_t e)
{
    slong anc, k;
    ulong v;

    v = fbnd_reduce(&e);
    anc = x->exp - x->size;
    k = e.a - anc;          /* e = v B^k units of the bottom limb */

    if (k > 0)
    {
        if (k > x->size - 1)
        {
            /* the bound swamps the mantissa: 0 +/- (|x| + e) with
               |x| < B^exp, at the anchor B^exp or above */
            fbnd_t t = fbnd_add(fbnd(0, v, e.a), fbnd(0, 1, x->exp));
            v = fbnd_reduce(&t);
            x->size = 0;
            x->negative = 0;
            x->exp = t.a;
            x->err = v;
            return;
        }
        /* drop k limbs, worth one unit */
        flint_mpn_copyi(x->d, x->d + k, x->size - k);
        x->size -= k;
        v++;
        if (v == 0)
        {
            /* the count was B - 1: one more limb, two units */
            flint_mpn_copyi(x->d, x->d + 1, x->size - 1);
            x->size -= 1;
            v = 2;
        }
    }
    else if (k < 0)
    {
        /* pad -k zero limbs below the mantissa */
        fball_fit(x, x->size - k);
        memmove(x->d - k, x->d, x->size * sizeof(ulong));
        flint_mpn_zero(x->d, -k);
        x->size -= k;
    }

    x->err = v;

    if (x->size > 0 && x->d[x->size - 1] == 0)
        _fball_norm(x);
}

FLINT_FORCE_INLINE void
_fball_apply_bnd(fball_t x, fbnd_t e)
{
    /* the common cases: an exact result, or a count below B already
       at the scale of the mantissa's bottom limb */
    if (e.hi == 0)
    {
        if (e.lo == 0)
        {
            x->err = 0;
            return;
        }
        if (x->size > 0 && e.a == x->exp - x->size)
        {
            x->err = e.lo;
            return;
        }
    }
    _fball_apply_bnd_slow(x, e);
}

static void
_fball_apply_err(fball_t x, fberr_t e)
{
    _fball_apply_bnd(x, fbnd_of_fberr(e));
}

/* set x to the ball 0 +/- e */
static void
_fball_zero_bnd(fball_t x, fbnd_t e)
{
    ulong v = fbnd_reduce(&e);
    x->size = 0;
    x->negative = 0;
    x->exp = e.a;
    x->err = v;
}

static void
_fball_zero_err(fball_t x, fberr_t e)
{
    _fball_zero_bnd(x, fbnd_of_fberr(e));
}

void
fball_init(fball_t x)
{
    x->d = NULL;
    x->alloc = 0;
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0;
}

void
fball_clear(fball_t x)
{
    flint_free(x->d);
}

void
_fball_grow(fball_t x, slong k)
{
    slong newalloc = FLINT_MAX(k, x->alloc + x->alloc / 2);
    x->d = flint_realloc(x->d, newalloc * sizeof(ulong));
    x->alloc = newalloc;
}

void
fball_zero(fball_t x)
{
    x->size = 0;
    x->negative = 0;
    x->exp = 0;
    x->err = 0;
}

int
fball_is_zero_exact(const fball_t x)
{
    return x->size == 0 && x->err == 0;
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
    x->err = 0;
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
}

void
fball_swap(fball_t x, fball_t y)
{
    FLINT_SWAP(fball_struct, *x, *y);
}

void
fball_add_error_ulps(fball_t x, double e)
{
    _fball_apply_bnd(x, fbnd_add(fbnd_of(x),
        fbnd_of_fberr(fberr(e, x->exp - x->size))));
}

void
fball_add_error(fball_t x, double v, slong anc)
{
    _fball_apply_bnd(x, fbnd_add(fbnd_of(x), fbnd_of_fberr(fberr(v, anc))));
}

void
fball_add_error_2exp_rel(fball_t x, slong e2)
{
    /* |value| < B^exp; add 2^e2 B^exp =
       2^(e2 mod FLINT_BITS) B^(exp + e2/FLINT_BITS), flooring the
       limb division so the bit remainder is >= 0 */
    slong q = e2 >> FBALL_LGB;
    int r = (int) (e2 - (q << FBALL_LGB));
    _fball_apply_bnd(x, fbnd_add(fbnd_of(x), fbnd(0, UWORD(1) << r, x->exp + q)));
}

/* ACCURATE LIMBS.  With errors normalized below 2^69 ulps, an
   inexact mantissa is accurate to all but its ~2 lowest limbs, so
   the mantissa length itself tracks accuracy; the target length of
   an output whose relative accuracy is limited by operand x is
   size + 1 when x is inexact (one guard limb; anything longer only
   computes noise) and unbounded when x is exact. */
/* THE NOISE FLOOR of an inexact ball is its bottom limb: the radius
   is a count of at least one ulp of that limb, so no ball holds limbs
   below it, and any limb an operation computes below the floor of its
   result is dropped again when the radius is installed.  (The
   double-anchored radius of an earlier design could sit far below one
   ulp, and operations kept a pad of six limbs under the floor; in the
   limb design the pad only cost a pass over the operand -- a sum whose
   other operand reached below an inexact one moved it down and back
   up -- and six limbs of every product.)

   A sum therefore takes nothing below an inexact operand's floor (the
   discarded tail of the other costs one unit there), and products and
   quotients of inexact operands keep size + 2 limbs: the truncation
   below the second guard limb is then 1/B of the inherited radius. */
static slong
_fball_noise_floor(const fball_t x)
{
    return x->exp - x->size;
}

#define FBALL_ACC_GUARD 2

static slong
_fball_acc(const fball_t x, slong n)
{
    if (x->err == 0)
        return WORD_MAX;
    return x->size + FBALL_ACC_GUARD;
}

/* the windowed middle product beats the full product from this many
   kept limbs on (measured on x86-64: at 32 limbs the window costs 1.3
   times the full product, at 64 the two are even, at 256 the window
   is 30% faster) */
#ifndef FBALL_MULMID_CUTOFF
#define FBALL_MULMID_CUTOFF 64
#endif

/* the radius of a product: the cross terms |a| errb + |b| erra +
   erra errb, valid also for zero mantissas, magnitudes through the
   leading bits (see _fball_lead: bounding through the top limb alone
   overstates by a factor up to 2 for a top limb of 1, and that
   pessimism compounds across a chain of multiplications) */
static fbnd_t
_fball_mul_bnd(const fball_t a, const fball_t b)
{
    fbnd_t e = fbnd_zero;

    if (b->err != 0)
    {
        if (a->size == 0)
            e = fbnd_add(e, fbnd(0, b->err, a->exp + b->exp - b->size));
        else
            e = fbnd_add(e, fbnd_mul_mag(b->err, a,
                    (a->exp - 1) + b->exp - b->size));
    }
    if (a->err != 0)
    {
        if (b->size == 0)
            e = fbnd_add(e, fbnd(0, a->err, b->exp + a->exp - a->size));
        else
            e = fbnd_add(e, fbnd_mul_mag(a->err, b,
                    (b->exp - 1) + a->exp - a->size));
        if (b->err != 0)
        {
            fbnd_t t;
            umul_ppmm(t.hi, t.lo, a->err, b->err);
            t.a = a->exp - a->size + b->exp - b->size;
            e = fbnd_add(e, t);
        }
    }
    return e;
}

/* res = a * b truncated to about n limbs.  Full products are used
   when at most a few limbs longer than the kept part (the low limbs
   then round down into <= 1 ulp) and below FBALL_MULMID_CUTOFF kept
   limbs, otherwise flint_mpn_mulmid computes just the kept window
   plus one guard limb below, whose one-sided deficit of up to
   min(an, bn) + 2 ulps one limb above the window bottom is absorbed
   into the radius (and immediately truncated away by the error
   normalization). */
void
fball_mul(fball_t res, const fball_t a, const fball_t b, slong n)
{
    slong sa = a->size, sb = b->size, ps, keep;
    slong rexp;
    int rneg, exact = (a->err == 0 && b->err == 0);
    fbnd_t e;
    nn_srcptr ad, bd;
    TMP_INIT;

    e = _fball_mul_bnd(a, b);

    if (sa == 0 || sb == 0)
    {
        _fball_zero_bnd(res, e);
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
        fball_fit(res, sz - lo);
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
        keep = FLINT_MIN(keep, _fball_acc(a, n));
        keep = FLINT_MIN(keep, _fball_acc(b, n));
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
        e = fbnd_add(e, fbnd_mul_mag(1, b,
                (a->exp - (keep + 2)) + (b->exp - 1)));
        ad += sa - (keep + 2);
        sa = keep + 2;
    }
    if (sb > keep + 2)
    {
        e = fbnd_add(e, fbnd_mul_mag(1, a,
                (b->exp - (keep + 2)) + (a->exp - 1)));
        bd += sb - (keep + 2);
        sb = keep + 2;
    }
    ps = sa + sb;

    if (ps <= keep + 5 || ps <= keep + keep / 16
        || keep < FBALL_MULMID_CUTOFF || (ad == bd && sa == sb))
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
            fball_fit(res, ps);
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
                e = fbnd_add(e, fbnd(0, 1, rexp - ps + drop));
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
        res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_bnd(res, e);
        TMP_END;
    }
    else
    {
        /* windowed middle product: window [ps - keep - 1, ps), one
           limb more when the top limb of the product is zero (short of
           a carry from below) */
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
        zlo = ps - keep - 1;
        zn = keep + 1;

        if (aliased)
            t = TMP_ALLOC(zn * sizeof(ulong));
        else
        {
            fball_fit(res, zn);
            t = res->d;
        }

        flint_mpn_mulmid(t, ad, sa, bd, sb, zlo, ps);

        /* deficit (below the true window) plus the discarded low
           limbs: < (min + 2) B + 1 ulps of the window bottom */
        e = fbnd_add(e, fbnd(0, FLINT_MIN(FLINT_MIN(sa, sb), zlo) + 2,
            rexp - ps + zlo + 1));
        e = fbnd_add(e, fbnd(0, 1, rexp - ps + zlo));

        if (aliased)
        {
            fball_fit(res, zn);
            flint_mpn_copyi(res->d, t, zn);
        }
        res->size = zn;
        res->exp = rexp;
        res->negative = rneg;
        res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_bnd(res, e);
        TMP_END;
    }
}

/* The window of one part of a complex operand at the common bottom
   limb bot of the operand: the part's limbs above bot (a pointer into
   its mantissa when its bottom lies at or below bot, else the mantissa
   over zero padding written into buf), at least one limb */
static nn_srcptr
_fball_cwin(const fball_t x, slong bot, nn_ptr buf, slong * len)
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
_fball_set_product(fball_t res, nn_srcptr z, slong zl, slong bot,
    slong cut, fbnd_t e)
{
    slong len = FLINT_ABS(zl), drop;

    while (len > 0 && z[len - 1] == 0)
        len--;
    drop = FLINT_MAX(0, cut - bot);
    if (drop >= len)
    {
        if (len > 0)
            e = fbnd_add(e, fbnd(0, 1, bot + len));
        if (fbnd_is_zero(e))
        {
            res->size = 0;
            res->negative = 0;
            res->exp = cut;
            res->err = 0;
        }
        else
            _fball_zero_bnd(res, e);
        return;
    }
    if (drop > 0)
    {
        slong i;
        int inexact = 0;
        for (i = 0; i < drop && !inexact; i++)
            inexact = (z[i] != 0);
        if (inexact)
            e = fbnd_add(e, fbnd(0, 1, bot + drop));
    }
    fball_fit(res, len - drop);
    flint_mpn_copyi(res->d, z + drop, len - drop);
    res->size = len - drop;
    res->exp = bot + len;
    res->negative = (zl < 0);
    res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_bnd(res, e);
}

/* |x| < B^mexp(x), for a zero mantissa through the radius */
static slong
_fball_mexp(const fball_t x)
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
fball_mul_complex(fball_t rr, fball_t ri, const fball_t ar,
    const fball_t ai, const fball_t br, const fball_t bi, slong n)
{
    slong topa, topb, T, fl, keep, bota, botb, cut, la[2], lb[2], zl;
    slong ea_r, ea_i, eb_r, eb_i, zr_len, zi_len;
    nn_srcptr wa[2], wb[2];
    nn_ptr buf, zr, zi;
    fbnd_t er, ei;
    TMP_INIT;

    ea_r = _fball_mexp(ar); ea_i = _fball_mexp(ai);
    eb_r = _fball_mexp(br); eb_i = _fball_mexp(bi);

    if ((ar->size == 0 && ai->size == 0) || (br->size == 0 && bi->size == 0))
    {
        /* a zero operand (possibly with a radius): the plain formulas */
        fball_t t1, t2;
        fball_init(t1);
        fball_init(t2);
        fball_mul(t1, ar, br, n);
        fball_mul(t2, ai, bi, n);
        fball_sub(t1, t1, t2, n);
        fball_mul(t2, ar, bi, n);
        fball_mul(ri, ai, br, n);
        fball_add(ri, ri, t2, n);
        fball_swap(rr, t1);
        fball_clear(t1);
        fball_clear(t2);
        return;
    }

    /* the frames: the tops of the operands and of the result, the
       noise floor of the result (the term p q with p inexact is noisy
       from mexp(q) + floor(p) up), and the kept limbs */
    topa = FLINT_MAX(ar->size ? ar->exp : WORD_MIN, ai->size ? ai->exp : WORD_MIN);
    topb = FLINT_MAX(br->size ? br->exp : WORD_MIN, bi->size ? bi->exp : WORD_MIN);
    T = topa + topb;
    fl = WORD_MIN;
    if (ar->err) fl = FLINT_MAX(fl, FLINT_MAX(eb_r, eb_i) + _fball_noise_floor(ar));
    if (ai->err) fl = FLINT_MAX(fl, FLINT_MAX(eb_r, eb_i) + _fball_noise_floor(ai));
    if (br->err) fl = FLINT_MAX(fl, FLINT_MAX(ea_r, ea_i) + _fball_noise_floor(br));
    if (bi->err) fl = FLINT_MAX(fl, FLINT_MAX(ea_r, ea_i) + _fball_noise_floor(bi));
    keep = n;
    if (fl != WORD_MIN)
        keep = FLINT_MIN(keep, T - fl + FBALL_ACC_GUARD);
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
    er = fbnd_add(_fball_mul_bnd(ar, br), _fball_mul_bnd(ai, bi));
    ei = fbnd_add(_fball_mul_bnd(ar, bi), _fball_mul_bnd(ai, br));
    if (ar->size && ar->exp - ar->size < bota)
    {
        if (br->size) er = fbnd_add(er, fbnd_mul_mag(1, br, bota + br->exp - 1));
        if (bi->size) ei = fbnd_add(ei, fbnd_mul_mag(1, bi, bota + bi->exp - 1));
    }
    if (ai->size && ai->exp - ai->size < bota)
    {
        if (bi->size) er = fbnd_add(er, fbnd_mul_mag(1, bi, bota + bi->exp - 1));
        if (br->size) ei = fbnd_add(ei, fbnd_mul_mag(1, br, bota + br->exp - 1));
    }
    if (br->size && br->exp - br->size < botb)
    {
        if (ar->size) er = fbnd_add(er, fbnd_mul_mag(1, ar, botb + ar->exp - 1));
        if (ai->size) ei = fbnd_add(ei, fbnd_mul_mag(1, ai, botb + ai->exp - 1));
    }
    if (bi->size && bi->exp - bi->size < botb)
    {
        if (ai->size) er = fbnd_add(er, fbnd_mul_mag(1, ai, botb + ai->exp - 1));
        if (ar->size) ei = fbnd_add(ei, fbnd_mul_mag(1, ar, botb + ar->exp - 1));
    }

    TMP_START;
    if (keep >= FBALL_MULMID_CUTOFF
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
            w = _fball_cwin(ar, bota, zr, &l);
            flint_mpn_zero(pa[0], ka); flint_mpn_copyi(pa[0] + ka, w, l); flint_mpn_zero(pa[0] + ka + l, nn - ka - l);
            w = _fball_cwin(ai, bota, zr, &l);
            flint_mpn_zero(pa[1], ka); flint_mpn_copyi(pa[1] + ka, w, l); flint_mpn_zero(pa[1] + ka + l, nn - ka - l);
            w = _fball_cwin(br, botb, zr, &l);
            flint_mpn_zero(pb[0], kb); flint_mpn_copyi(pb[0] + kb, w, l); flint_mpn_zero(pb[0] + kb + l, nn - kb - l);
            w = _fball_cwin(bi, botb, zr, &l);
            flint_mpn_zero(pb[1], kb); flint_mpn_copyi(pb[1] + kb, w, l); flint_mpn_zero(pb[1] + kb + l, nn - kb - l);
        }

        flint_mpn_mulhigh_n_complex(zr, &sr, zi, &si,
            pa[0], ar->negative, pa[1], ai->negative,
            pb[0], br->negative, pb[1], bi->negative, nn);

        er = fbnd_add(er, fbnd(0, 5, scale));
        ei = fbnd_add(ei, fbnd(0, 5, scale));
        zl = nn + 1;
        _fball_set_product(rr, zr, sr ? -zl : zl, scale, cut, er);
        _fball_set_product(ri, zi, si ? -zl : zl, scale, cut, ei);
        TMP_END;
        return;
    }

    zl = FLINT_MAX(la[0], la[1]) + FLINT_MAX(lb[0], lb[1]) + 1;
    buf = TMP_ALLOC((la[0] + la[1] + lb[0] + lb[1] + 2 * zl) * sizeof(ulong));
    wa[0] = _fball_cwin(ar, bota, buf, &la[0]);
    wa[1] = _fball_cwin(ai, bota, buf + la[0], &la[1]);
    wb[0] = _fball_cwin(br, botb, buf + la[0] + la[1], &lb[0]);
    wb[1] = _fball_cwin(bi, botb, buf + la[0] + la[1] + lb[0], &lb[1]);
    zr = buf + la[0] + la[1] + lb[0] + lb[1];
    zi = zr + zl;

    flint_mpn_mul_complex(zr, &zr_len, zi, &zi_len,
        wa[0], la[0], ar->negative, wa[1], la[1], ai->negative,
        wb[0], lb[0], br->negative, wb[1], lb[1], bi->negative);

    _fball_set_product(rr, zr, zr_len, bota + botb, cut, er);
    _fball_set_product(ri, zi, zi_len, bota + botb, cut, ei);
    TMP_END;
}

void
fball_mul_ui(fball_t res, const fball_t a, ulong c, slong n)
{
    fbnd_t e;
    nn_ptr t;
    ulong cy;
    slong sa = a->size, keep, drop;

    if (c == 0)
    {
        fball_zero(res);
        return;
    }

    e = fbnd_zero;
    if (a->err != 0)
    {
        umul_ppmm(e.hi, e.lo, a->err, c);
        e.a = (sa == 0) ? a->exp : a->exp - a->size;
    }

    if (sa == 0)
    {
        _fball_zero_bnd(res, e);
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
                e = fbnd_add(e, fbnd(0, 1, anc + drop));
            flint_mpn_copyi(t, t + drop, keep);
        }

        res->size = keep;
        res->exp = anc + ts;
        res->negative = a->negative;
        res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
        _fball_norm(res);
        _fball_apply_bnd(res, e);
    }
}

/* v B / c rounded up as a count at the anchor below v's: v exact in
   a double below 2^53 (c within 2^-53, the quotient within 2^-52, so
   a factor 1 + 2^-50 and one unit up cover the rounding), else by an
   integer division of the rare wide radius */
FLINT_FORCE_INLINE fbnd_t
_fbnd_div_ui_B(ulong v, ulong c, slong a)
{
    fbnd_t e;
    if (v < (UWORD(1) << FLINT_MIN(53, FLINT_BITS - 1)))
    {
        double d = (double) v / (double) c * (1.0 + 0x1p-50), f;
        e.hi = (ulong) d;
        f = (d - (double) e.hi) * FBALL_D_B;
        e.lo = (f >= FBALL_D_B - 2.0) ? (ulong) (FBALL_D_B - 2.0) : (ulong) f;
        add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(2));
    }
    else
    {
        ulong q, r;
        q = v / c;
        r = v - q * c;
        udiv_qrnnd(e.lo, r, r, UWORD(0), c);
        e.hi = q;
        add_ssaaaa(e.hi, e.lo, e.hi, e.lo, UWORD(0), UWORD(1));
    }
    e.a = a;
    return e;
}

/* res = a / c truncated to about n limbs, c >= 1: the quotient of the
   top min(n + 2, size) limbs of the mantissa over one limb of zeros
   (two for an exact dividend) by mpn_divrem_1, so that the result keeps
   the quotient's own precision: a quotient re-anchored at the
   dividend's bottom limb with its truncation error there loses
   log2(c) bits relative to its magnitude at every division, which a
   chain of divisions (a common denominator k^m) compounds to nothing.
   The radius err / c and the dropped low limbs of the dividend are
   bounded one limb below the dividend's anchor, where the truncation
   of the quotient is one unit; for an exact dividend the quotient
   runs two limbs below, its truncation one unit there. */
void
fball_div_ui(fball_t res, const fball_t a, ulong c, slong n)
{
    fbnd_t e;
    nn_ptr t;
    slong sa = a->size, keep, drop, anc, qxn;

    if (c == 0)
        flint_throw(FLINT_ERROR, "fball_div_ui: division by zero\n");

    if (sa == 0)
    {
        e = (a->err != 0) ? _fbnd_div_ui_B(a->err, c, a->exp - 1) : fbnd_zero;
        _fball_zero_bnd(res, e);
        return;
    }

    keep = FLINT_MIN(n + 2, sa);
    drop = sa - keep;
    anc = a->exp - a->size + drop;  /* ulp of the kept part */

    if (a->err == 0 && drop == 0)
    {
        qxn = 2;
        e = fbnd(0, 1, anc - 2);
    }
    else
    {
        /* the radius and the dropped limbs of a (below B^anc, so below
           B^anc / c after division) at B^(anc - 1), the quotient's
           truncation one unit there */
        qxn = 1;
        e = _fbnd_div_ui_B(a->err + (drop > 0), c, anc - 1);
    }

    fball_fit(res, keep + qxn);
    t = res->d;
    if (res != a)
        mpn_divrem_1(t, qxn, a->d + drop, keep, c);
    else
    {
        memmove(t + qxn, t + drop, keep * sizeof(ulong));
        t[0] = 0;
        if (qxn == 2)
            t[1] = 0;
        mpn_divrem_1(t, 0, t, keep + qxn, c);
    }

    res->size = keep + qxn;
    res->exp = anc + keep;
    res->negative = a->negative;
    res->err = 1;       /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_bnd(res, e);
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
   the operands enter as in fball_mul.  With the assumption violated
   the result is wrong (the window wraps), so it must be rigorous: a
   ball computed by an ordinary fball_mul and fball_sub would have been
   correct at any magnitude, at the cost of a high product of b c that
   is known to cancel. */
void
fball_submul_bounded(fball_t res, const fball_t a, const fball_t b,
    const fball_t c, slong E, slong n)
{
    slong sb = b->size, sc = c->size, sa = a->size;
    slong wlo = E - n - 2, wn = n + 3, base, zlo, zhi, i, lo, hi;
    nn_srcptr bd, cd;
    nn_ptr t;
    fbnd_t e;
    int neg;
    TMP_INIT;

    /* the radii: err_a + |b| err_c + |c| err_b + err_b err_c, through
       the top limbs as in fball_mul */
    e = fbnd_of(a);
    if (c->err != 0)
    {
        if (sb == 0)
            e = fbnd_add(e, fbnd(0, c->err, b->exp + c->exp - c->size));
        else
            e = fbnd_add(e, fbnd_mul_mag(c->err, b,
                    (b->exp - 1) + c->exp - c->size));
    }
    if (b->err != 0)
    {
        if (sc == 0)
            e = fbnd_add(e, fbnd(0, b->err, c->exp + b->exp - b->size));
        else
            e = fbnd_add(e, fbnd_mul_mag(b->err, c,
                    (c->exp - 1) + b->exp - b->size));
        if (c->err != 0)
        {
            fbnd_t t2;
            umul_ppmm(t2.hi, t2.lo, b->err, c->err);
            t2.a = b->exp - b->size + c->exp - c->size;
            e = fbnd_add(e, t2);
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
    fball_fit(res, n);
    flint_mpn_copyi(res->d, t + 2, n);
    for (i = n; i > 0 && res->d[i - 1] == 0; i--)
        ;
    res->size = i;
    res->exp = wlo + 2 + i;
    res->negative = neg;
    e = fbnd_add(e, fbnd(0, 3, wlo + 2));
    res->err = 1;
    if (res->size == 0)
        _fball_zero_bnd(res, e);
    else
        _fball_apply_bnd(res, e);

    TMP_END;
}

/* zero a few limbs without a memset call */
FLINT_FORCE_INLINE void
_fball_zerow(nn_ptr p, slong k)
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
_fball_addw(nn_ptr t, slong tn, nn_srcptr b, slong bn)
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
_fball_subw(nn_ptr t, slong tn, nn_srcptr b, slong bn)
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
_fball_add_short(fball_t res, const fball_t a, int aneg, const fball_t b,
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
    fball_fit(res, k - i);
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
_fball_mag_ge(const fball_t a, const fball_t b)
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
_fball_add(fball_t res, const fball_t a, int aneg, const fball_t b,
    int bneg, slong n)
{
    slong top, bot, span, oa, ob;
    fbnd_t e;
    nn_ptr t;
    int rneg, rev = 0, trunc_a = 0, trunc_b = 0;

    if (_fball_add_short(res, a, aneg, b, bneg, n))
        return;

    if ((a->err | b->err) == 0)
        e = fbnd_zero;
    else
        e = fbnd_add(fbnd_of(a), fbnd_of(b));

    if (a->size == 0 || b->size == 0)
    {
        if (a->size == 0 && b->size == 0)
        {
            _fball_zero_bnd(res, e);
            return;
        }
        /* e already contains both radii; apply_bnd overwrites err */
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
        _fball_apply_bnd(res, e);
        return;
    }

    if (a == b)
    {
        if (aneg == bneg)
        {
            fball_mul_ui(res, a, 2, n);
            res->negative = aneg;
        }
        else
            _fball_zero_bnd(res, e);
        return;
    }

    if (res == b)
    {
        const fball_struct * s = a; a = b; b = s;
        rneg = aneg; aneg = bneg; bneg = rneg;
    }

    if (aneg != bneg && !_fball_mag_ge(a, b))
    {
        if (res != a)
        {
            const fball_struct * s = a; a = b; b = s;
            rneg = aneg; aneg = bneg; bneg = rneg;
        }
        else
            rev = 1;
    }

    top = FLINT_MAX(a->exp, b->exp) + 1;
    bot = FLINT_MIN(a->exp - a->size, b->exp - b->size);
    bot = FLINT_MAX(bot, top - (n + 2));
    /* nothing below an inexact operand's noise floor (see
       _fball_noise_floor): the other operand's tail there is dropped
       for one unit of the floor, and an in-place operand never moves */
    if (a->err != 0)
        bot = FLINT_MAX(bot, _fball_noise_floor(a));
    if (b->err != 0)
        bot = FLINT_MAX(bot, _fball_noise_floor(b));

    span = top - bot;
    oa = (a->exp - a->size) - bot;

    fball_fit(res, span);
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
            _fball_zerow(t, oa);
        _fball_zerow(t + oa + a->size, span - oa - a->size);
    }
    else if (a->size + oa > 0)
    {
        /* low limbs of a truncated; copying downward is safe in
           place (ascending copy, dst below src) */
        flint_mpn_copyi(t, (res == a ? t : a->d) - oa, a->size + oa);
        _fball_zerow(t + a->size + oa, span - a->size - oa);
        trunc_a = 1;
    }
    else
    {
        trunc_a = 1;
        _fball_zerow(t, span);
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
                ulong cy = _fball_addw(t + ob, span - ob, bd, bs);
                FLINT_ASSERT(cy == 0);  /* top slot reserved */
                (void) cy;
            }
            else if (!rev)
            {
                ulong bw = _fball_subw(t + ob, span - ob, bd, bs);
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
        e = fbnd_add(e, fbnd(0, 1, bot));
    if (trunc_b)
        e = fbnd_add(e, fbnd(0, 1, bot));

    /* the reserved carry slot is usually empty */
    if (t[span - 1] == 0)
    {
        span--;
        top--;
    }
    res->size = span;
    res->exp = top;
    res->negative = rneg;
    res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
    _fball_norm(res);
    if (res->size == 0)
        _fball_zero_bnd(res, e);
    else
        _fball_apply_bnd(res, e);
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

/* res = x + (-1)^sub y c.  In place (res == x) whenever the product
   lands inside x's window -- y c below x's top limb, x no longer than
   the precision asks, for a difference |y c| < B^(x.exp - 1) <= |x|
   so that the sign cannot flip, and y's bottom at or above x's, or
   anywhere below it when x is inexact (its noise floor is its bottom
   limb, see _fball_noise_floor) -- by mpn_addmul_1 resp.
   mpn_submul_1 over y's limbs from x's bottom up and a carry or
   borrow that stops when it runs out: O(|y|), whatever the size of x.
   The dropped tail of y, below one unit of x's bottom limb, is below
   c units there after the multiplication: the high word of c times
   its top limb enters as a carry, and the rest is under 2 units.  The
   radius is x's plus c times y's (exact products of words), added as
   in fball_add.  Otherwise (and when res aliases y) y c is formed by
   fball_mul_ui and added by fball_add. */
static void
_fball_addmul_ui(fball_t res, const fball_t x, const fball_t y, ulong c,
    slong n, int sub)
{
    slong xbot, ybot, off, xs, ys, drop;
    int ysgn;
    fbnd_t e;

    if (c == 0 || fball_is_zero_exact(y))
    {
        fball_set(res, x);
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
            fball_set(res, x);
        d = res->d;

        e = fbnd_of(res);
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
            e = fbnd_add(e, fbnd(0, 2, xbot));
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
                fball_fit(res, xs + 1);
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
            fbnd_t t;
            umul_ppmm(t.hi, t.lo, y->err, c);
            t.a = ybot;
            e = fbnd_add(e, t);
        }
        res->err = !fbnd_is_zero(e);    /* norm must know inexactness */
        _fball_norm(res);
        if (res->size == 0)
            _fball_zero_bnd(res, e);
        else
            _fball_apply_bnd(res, e);
        return;
    }

    {
        fball_t t;
        fball_init(t);
        fball_mul_ui(t, y, c, n);
        if (sub)
            fball_sub(res, x, t, n);
        else
            fball_add(res, x, t, n);
        fball_clear(t);
    }
}

void
fball_addmul_ui(fball_t res, const fball_t x, const fball_t y, ulong c,
    slong n)
{
    _fball_addmul_ui(res, x, y, c, n, 0);
}

void
fball_submul_ui(fball_t res, const fball_t x, const fball_t y, ulong c,
    slong n)
{
    _fball_addmul_ui(res, x, y, c, n, 1);
}

/* rigorous OVERestimate of the relative radius err / |x| of a ball
   with nonzero mantissa: |x| >= top B^(exp - 1), so
   rel <= (err / top) B^(1 - size); the ldexp window is clamped on
   the small side to a value exceeding anything the invariant
   err < B allows there, keeping the estimate one-sided */
static double
_fball_rel_bound(const fball_t x)
{
    slong k = 1 - x->size;

    if (x->err == 0)
        return 0.0;
    if (x->size == 0)
        return HUGE_VAL;
    if (FLINT_BITS * k > 900)
        return HUGE_VAL;
    if (FLINT_BITS * k < -900)
        return 0x1p-700;    /* true bound < 2^(FLINT_BITS - 900) */
    /* |x| >= top B^(exp - 1): through the top limb, as the operation
       bounds do (the blanket B^(exp - 1) overstates by up to a limb) */
    return ldexp((double) x->err / _fball_mag_lo(x),
        (int) (FLINT_BITS * k)) * FBALL_EPS;
}

/* res = a / b.  From FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF limbs (of
   both the divisor and the quotient) the mantissas are divided as
   fractions by fixed_div_newton, which takes the divisor as it is (no
   normalization shifts; one extra guard limb absorbs a small top
   limb); below, flint_mpn_divapprox_fraction picks the register,
   schoolbook, blockwise or divide-and-conquer division and needs
   neither zero-extended nor normalized operands.  Operand errors
   enter as err_a / |b| + |a| err_b / |b|^2, with |a| and |b| bounded
   through their top limbs and the -err_b correction of the
   denominator covered by a 1.01 factor, which is rigorous because the
   relative radius of b is CHECKED to be below 2^-30 (a divisor ball
   this wide is a usage error: the mantissa would be pure noise). */
#define FBALL_DIV_NEWTON_CUTOFF FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF

void
fball_div(fball_t res, const fball_t a, const fball_t b, slong n)
{
    slong sa = a->size, sb = b->size, nd, qn;
    int rneg;
    double mb;
    fberr_t e;
    TMP_INIT;

    if (sb == 0)
        flint_throw(FLINT_ERROR, "fball_div: division by zero ball\n");
    if (_fball_rel_bound(b) > 0x1p-30)
        flint_throw(FLINT_ERROR,
            "fball_div: divisor relative radius above 2^-30\n");

    /* magnitudes through the top limbs, as in fball_mul:
       |a| < ma B^(a.exp - 1), |b| >= mb B^(b.exp - 1) with
       ma = a_top + 1 rounded up and mb = b_top rounded down; the
       blanket bounds |a| < B^a.exp, |b| >= B^(b.exp - 1) overstate
       the quotient by up to a limb each, and |a| / |b|^2 by up to
       three */
    e = fberr(0.0, 0);
    mb = _fball_mag_lo(b);
    {
        if (a->err != 0)
            e = fberr_add(e, fberr((double) a->err * 1.01 / mb,
                    (a->exp - sa) + 1 - b->exp));
        if (b->err != 0)
        {
            double ma = (sa == 0) ? 1.0 : _fball_mag_hi(a);
            slong ea = (sa == 0) ? a->exp + 1 : a->exp;
            e = fberr_add(e, fberr((double) b->err * 1.01 * (ma / (mb * mb)),
                    ea - 1 + (b->exp - sb) + 2
                        - 2 * b->exp));
        }
    }

    if (sa == 0)
    {
        _fball_zero_err(res, e);
        return;
    }

    nd = FLINT_MIN(n, FLINT_MIN(_fball_acc(a, n), _fball_acc(b, n)));
    nd = FLINT_MAX(nd, 2);
    qn = nd + 2;

    rneg = a->negative ^ b->negative;

    TMP_START;
    {
        int aliased = (res == a || res == b);
        slong rn, E;
        nn_ptr q;

        if (sb >= FBALL_DIV_NEWTON_CUTOFF && qn >= FBALL_DIV_NEWTON_CUTOFF)
        {
            /* Karp-Markstein on the mantissas read as fractions,
               (a_d B^-sa) / (b_d B^-sb): fixed_div_newton does not
               need a normalized divisor, only a nonzero top limb, so
               den = b_d B^-sb >= b_top / B and the error
               4 B^-(nd+1) / den is at most 4 B / b_top ulps of the
               nd + 1 fraction limbs -- one more than otherwise needed,
               which covers the unnormalized divisor */
            rn = nd + 3;
            q = aliased ? TMP_ALLOC(rn * sizeof(ulong))
                : (fball_fit(res, rn), res->d);
            fixed_div_newton(q, a->d, sa, b->d, sb, nd + 1);
            E = a->exp - b->exp + 2 - rn;
            e = fberr_add(e, fberr(4.0 / mb, E + 1));
        }
        else
        {
            /* approximate division by the mantissa of b (truncated
               to qn + 2 limbs when longer): Q = floor(a_d B^f / b_d)
               or one more, with f chosen for qn + 1 quotient limbs
               (negative f truncates the numerator, exactly as a
               floor); |Q - true| < 1 ulp of B^E */
            nn_srcptr bd = b->d;
            slong sbt = sb, f;

            if (sb > qn + 2)
            {
                /* b' = b - delta, 0 <= delta < B^(b.exp - sbt):
                   a/b' - a/b <= |a| delta / b'^2
                   <= (ma / mb^2) B^(a.exp - b.exp - sbt + 1) */
                double ma = _fball_mag_hi(a);
                sbt = qn + 2;
                bd = b->d + (sb - sbt);
                e = fberr_add(e, fberr(1.01 * ma / (mb * mb),
                    a->exp - b->exp - sbt + 1));
            }

            rn = qn + 1;
            f = qn - sa + sbt;
            q = aliased ? TMP_ALLOC(rn * sizeof(ulong))
                : (fball_fit(res, rn), res->d);
            flint_mpn_divapprox_fraction(q, a->d, sa, bd, sbt, f);
            E = (a->exp - sa) - f - (b->exp - sbt);
            e = fberr_add(e, fberr(1.0, E));
        }

        if (aliased)
        {
            fball_fit(res, rn);
            flint_mpn_copyi(res->d, q, rn);
        }

        res->size = rn;
        res->exp = E + rn;
    }
    res->negative = rneg;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _fball_norm(res);
    if (res->size == 0)
        _fball_zero_err(res, e);
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
    res->err = 1;       /* norm must know inexactness */
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
    x->err = 0;
    _fball_norm(x);
}

double
fball_get_fixed(nn_ptr y, slong wn, const fball_t x)
{
    double e;
    slong sh;

    /* radius in output ulps: err B^(exp - size + wn), saturating.
       Both saturation directions must OVERestimate: far above, the
       bound is useless anyway (HUGE_VAL); far below, err < B
       at k <= -k0 limbs is at most 2^(FLINT_BITS (1 - k0)), so the
       clamp constant must exceed that, while the ldexp in the live
       window (|FLINT_BITS k| <= 900) cannot underflow */
    if (x->err == 0)
        e = 0.0;
    else
    {
        slong k = x->exp - x->size + wn;
        if (FLINT_BITS * k > 900)
            e = HUGE_VAL;
        else if (FLINT_BITS * k < -900)
            e = 0x1p-700;   /* true bound < 2^(64 - 900) */
        else
            e = ldexp((double) x->err, (int) (FLINT_BITS * k)) * FBALL_EPS;
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
   the pure power of two 2^s, s = the bit length of err -- must fit
   strictly inside the tail on both
   sides: 2^s <= L (the true value cannot borrow below the grid line)
   and L < B^t - 2^s (nor carry above the next), both checked by bit
   scans without materializing L + 2^s. */
int
fball_get_fixed_floor(nn_ptr y, slong n, const fball_t x)
{
    slong sh, t, s, i, hi;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(x->exp <= 0);
    FLINT_ASSERT(x->size == 0 || !x->negative);

    flint_mpn_zero(y, n);

    if (x->size == 0)
    {
        /* 0 +/- err: the floor is 0 iff the radius stays strictly
           below B^-n: err B^exp < B^-n */
        if (x->err == 0)
            return 1;
        return (slong) FLINT_BIT_COUNT(x->err) + FLINT_BITS * (x->exp + n) <= 0;
    }

    sh = x->exp - x->size + n;  /* x B^n = (d, size) B^sh */
    t = -sh;

    if (t <= 0)
    {
        /* mantissa entirely on or above the grid: exact values
           only (any radius straddles the grid line the value sits
           on) */
        if (x->err != 0)
            return 0;
        FLINT_ASSERT(x->size + sh <= n);
        flint_mpn_copyi(y + sh, x->d, x->size);
        return 1;
    }

    if (t < x->size)
        flint_mpn_copyi(y, x->d + t, x->size - t);
    /* else the whole value sits below the grid: floor 0, and the
       t-limb tail is the mantissa padded above with zero limbs */

    if (x->err == 0)
        return 1;

    s = (slong) FLINT_BIT_COUNT(x->err);

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
    if (x->size == 0)
        arf_zero(arb_midref(res));
    else
    {
        /* the mantissa is exact at size limbs; value d B^(exp - size)
           has arf exponent FLINT_BITS exp + (bits shifted out) */
        slong fix;
        _arf_set_round_mpn(arb_midref(res), &fix, x->d, x->size,
            x->negative, FLINT_BITS * x->size, ARF_RND_DOWN);
        fmpz_set_si(ARF_EXPREF(arb_midref(res)),
            FLINT_BITS * x->exp + fix);
    }

    if (x->err == 0)
        mag_zero(arb_radref(res));
    else
    {
        mag_set_ui(arb_radref(res), x->err);
        mag_mul_2exp_si(arb_radref(res), arb_radref(res),
            FLINT_BITS * (x->exp - x->size));
    }
}

/* the relative radius of x (nonzero mantissa) is below 2^e: through
   the top limb, err < 2^bits(err); -WORD_MAX/2 when exact */
slong
fball_rel_2exp(const fball_t x)
{
    if (x->err == 0)
        return -WORD_MAX / 2;
    return FLINT_BITS * (1 - x->size) + (slong) FLINT_BIT_COUNT(x->err)
        - (FLINT_BIT_COUNT(x->d[x->size - 1]) - 1) + 1;
}

/* x += [-2^e, 2^e] */
void
fball_add_error_2exp(fball_t x, slong e)
{
    slong q = e >> FBALL_LGB;
    fball_add_error(x, ldexp(1.0, (int) (e - q * FLINT_BITS)), q);
}

/* e with |x| < 2^e for the ball (value plus radius), -WORD_MAX/2 for
   an exact zero */
slong
fball_mag_2exp(const fball_t x)
{
    slong e;

    if (x->size == 0)
        e = -WORD_MAX / 2;
    else
        e = FLINT_BITS * (x->exp - 1)
            + FLINT_BIT_COUNT(x->d[x->size - 1]);

    if (x->err != 0)
        e = FLINT_MAX(e, FLINT_BITS * (x->exp - x->size)
            + (slong) FLINT_BIT_COUNT(x->err)) + 1;
    return e;
}

void
fball_print(const fball_t x)
{
    arb_t t;
    arb_init(t);
    fball_get_arb(t, x);
    flint_printf("[size %wd, exp %wd, err %wu] = ",
        x->size, x->exp, x->err);
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
        /* the radius scales exactly: err 2^r = err 2^(r - B) B at the
           anchor moved up by one limb when r != 0 */
        if (x->err != 0 && r != 0)
        {
            fbnd_t t = fbnd(x->err >> (FLINT_BITS - r), x->err << r, x->exp - 1);
            _fball_zero_bnd(x, t);
        }
        return;
    }

    x->exp += q;
    if (r != 0)
    {
        ulong cy, err = x->err;

        fball_fit(x, x->size + 1);
        cy = mpn_lshift(x->d, x->d, x->size, r);
        if (cy)
        {
            x->d[x->size] = cy;
            x->size++;
            x->exp++;
        }
        _fball_norm(x);
        if (err != 0)
            _fball_apply_bnd(x, fbnd(err >> (FLINT_BITS - r), err << r,
                x->exp - x->size));
        else
            x->err = 0;
    }
}

/* shared alignment for sqrt/rsqrt: express x = ahat * B^E with E
   even, ahat in [B^-2, 1) given as an an-limb fraction (aa, an);
   the input radius keeps its anchor: err ulps of B^(E - an) */
#define FBALL_SQRT_ALIGN                                        \
    slong E = x->exp, an = x->size, nd;                         \
    nn_ptr aa;                                                  \
    fberr_t e;                                                  \
    double ml;                                                  \
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

    /* operand error through the derivative: with |x| >= ml B^(E-1-q)
       (ml the leading bits of the mantissa, q = 1 iff the alignment
       padded an odd exponent with a zero top limb),
       |d(1/sqrt x)/dx| = 1/(2 x^(3/2)) <= ml^(-3/2) B^(3(1+q-E)/2)/2
       and the radius is err B^(E - an), so
       |Delta| <= 0.51 err ml^(-3/2) B^(3(1+q)/2) B^(-an-E/2). */
    e = fberr(0.0, 0);
    ml = _fball_mag_lo(x);
    if (x->err != 0)
    {
        int q = (int) (x->exp & 1);
        e = fberr((double) x->err * 0.51 * pow(ml, -1.5)
            * (q ? FBALL_D_B * FBALL_D_B * FBALL_D_B : FBALL_D_B * FBALL_D_SQRTB),
            -an - E / 2);
    }

    fball_fit(res, nd + 2);
    fixed_rsqrt_newton(res->d, aa, an, nd);

    /* Newton: <= 4 B^-nd / sqrt(ahat) with ahat = ml B^(-1-q), so
       <= 4 ml^(-1/2) B^((1+q)/2) B^-nd, scaled B^(-E/2) */
    {
        int q = (int) (x->exp & 1);
        e = fberr_add(e, fberr(4.0 * sqrt(1.0 / ml)
            * (q ? FBALL_D_B : FBALL_D_SQRTB), -nd - E / 2));
    }

    res->size = nd + 2;
    res->exp = 2 - E / 2;
    res->negative = 0;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_err(res, e);
    TMP_END;
}

/* fball_sqrt takes the integer square root flint_mpn_sqrtrem (no
   remainder) of a 2m-limb input while the m-limb result is below this,
   fixed_sqrt_newton above: the Newton code, which needs no exact
   remainder, overtakes at about half the input length where
   flint_mpn_sqrtrem switches to its own Newton code (measured 1.6-2.5x
   faster below, crossover at m ~ 1300-2000 on x86-64) */
#define FBALL_SQRT_NEWTON_CUTOFF (FLINT_MPN_SQRTREM_NEWTON_CUTOFF / 4)

/* res = sqrt(x); x > 0, n-limb target */
void
fball_sqrt(fball_t res, const fball_t x, slong n)
{
    FBALL_SQRT_ALIGN

    /* operand error through the derivative: with |x| >= ml B^(E-1-q)
       as above, |d(sqrt x)/dx| = 1/(2 sqrt x) <= ml^(-1/2) B^((1+q-E)/2)/2
       and the radius err B^(E - an):
       |Delta| <= 0.51 err ml^(-1/2) B^((1+q)/2) B^(E/2 - an) */
    e = fberr(0.0, 0);
    ml = _fball_mag_lo(x);
    if (x->err != 0)
    {
        int q = (int) (x->exp & 1);
        e = fberr((double) x->err * 0.51 * sqrt(1.0 / ml)
            * (q ? FBALL_D_B : FBALL_D_SQRTB), E / 2 - an);
    }

    if (nd + 2 < FBALL_SQRT_NEWTON_CUTOFF)
    {
        /* integer square root of X = ahat B^(2m) (the aligned
           mantissa zero-extended, or truncated, to 2m limbs):
           S = floor(sqrt(X)) has m limbs and sqrt(x) = S B^(E/2 - m)
           to within 1 ulp, plus < 1/2 ulp when X was truncated
           (sqrt(X) - sqrt(X - t) < t / (2 sqrt(X - t)) for t < 1) */
        slong m = nd + 2, xn = 2 * m, t = xn - an;
        nn_ptr X;

        X = TMP_ALLOC(xn * sizeof(ulong));
        if (t >= 0)
        {
            flint_mpn_zero(X, t);
            flint_mpn_copyi(X + t, aa, an);
        }
        else
        {
            flint_mpn_copyi(X, aa - t, xn);
        }
        while (X[xn - 1] == 0)      /* the alignment's zero top limb */
            xn--;

        fball_fit(res, m);
        flint_mpn_zero(res->d, m);
        flint_mpn_sqrtrem(res->d, NULL, X, xn);
        e = fberr_add(e, fberr(1.5, E / 2 - m));

        res->size = m;
        res->exp = E / 2;
    }
    else
    {
        /* the Newton error 4 B^-nd / sqrt(ahat) is absolute while
           sqrt(ahat) can be as small as B^-1 (an odd exponent aligned
           by a zero top limb) or B^-(1/2): one or two limbs more of
           the fraction keep the result accurate to nd limbs of its
           own */
        int q = (int) (x->exp & 1);
        slong nd2 = nd + 1 + q;

        fball_fit(res, nd2 + 2);
        fixed_sqrt_newton(res->d, aa, an, nd2);

        /* Newton: <= 4 B^-nd2 / sqrt(ahat) with ahat = ml B^(-1-q):
           <= 4 ml^(-1/2) B^((1+q)/2) B^-nd2, scaled B^(E/2) */
        e = fberr_add(e, fberr(4.0 * sqrt(1.0 / ml)
            * (q ? FBALL_D_B : FBALL_D_SQRTB), -nd2 + E / 2));

        res->size = nd2 + 2;
        res->exp = 2 + E / 2;
    }
    res->negative = 0;
    res->err = (e.v != 0.0);    /* norm must know inexactness */
    _fball_norm(res);
    _fball_apply_err(res, e);
    TMP_END;
}

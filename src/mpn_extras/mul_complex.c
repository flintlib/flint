/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpn_extras.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/*
    Complex multiplication of n-limb components with separate sign bits
    (0 = nonnegative), for numerical applications:

        zr + i zi = (ar + i ai) * (br + i bi)

    with all values given as (magnitude, sign bit).

    The full products are exact: zr and zi receive 2n + 1 limbs (the
    extra limb holds the carry of the additive combination). Below the
    fft threshold they use a Karatsuba-style scheme -- three real
    multiplications for a product, via t3 = (ar + ai)(br + bi), and two
    for a square, via zr = (ar + ai)(ar - ai) -- with the terms
    accumulated in a two's-complement accumulator so that the mixed signs
    need no case analysis. Above the threshold every variant (full, low,
    high and windowed products and squares) goes through one fft_small
    engine (_fft_small_mpn_mul_complex and relatives), which transforms
    each operand once and reuses it across the pointwise products:
    4 forward + 2 inverse transforms replace the 8 + 4 of four fused
    multiplications (2 + 2 replace 6 + 3 for squaring).

    The high variants receive n + 1 limbs, the limbs [n, 2n] of the exact
    result. They are approximate, as usual for high products: each real
    product is a mulhigh (or sqrhigh), whose window is exact or one unit
    low, and the discarded low halves can carry; above the fft threshold
    the truncated reconstruction perturbs its lowest limb by at most one
    unit. Summing the terms leaves an error of at most four units in the
    lowest returned limb (the worst observed over wide random testing is
    two) -- the caller is expected to track it, as ball arithmetic does.

    Karatsuba applies on the high path too. Its auxiliary product
    (ar + ai)(br + bi) has (n+1)-limb operands, so a plain mulhigh would
    start one limb above the window the other terms supply; giving each
    operand a low zero limb scales the product by 2^128, and the high
    n + 2 limbs of the scaled product are exactly the limbs [n, 2n+1] of
    the unscaled one -- the wanted window with a guard limb. That costs
    one mulhigh of size n + 2 in place of two of size n.

    */

slong flint_mpn_mul_complex_fft_cutoff = 230;
slong flint_mpn_sqr_complex_fft_cutoff = 270;

/* Defined here rather than in the fft_small module, which unsupported
   platforms exclude from the build entirely: the budget must be
   settable everywhere, since the drivers' fallback logic references it
   unconditionally. */
#if FLINT_BITS == 64
ulong flint_fft_small_max_transformed_ring_size = UWORD(4) << 30;
#else
ulong flint_fft_small_max_transformed_ring_size = UWORD(1) << 30;
#endif

#define MUL_COMPLEX_KARATSUBA_CUTOFF 8
#define SQR_COMPLEX_KARATSUBA_CUTOFF 16
#define MULHIGH_COMPLEX_KARATSUBA_CUTOFF 12
#define SQRHIGH_COMPLEX_KARATSUBA_CUTOFF 20

#if FLINT_HAVE_FFT_SMALL

/* the high variants promise exactly zlimbs limbs, zero padded, plus a
   sign; the shared fft path hands back a signed length */
static void
_high_pad(nn_ptr z, int * zsgn, slong len, mp_size_t zlimbs)
{
    mp_size_t rn = FLINT_ABS(len);

    if (rn < zlimbs)
        flint_mpn_zero(z + rn, zlimbs - rn);
    *zsgn = (len < 0);
}

#endif /* FLINT_HAVE_FFT_SMALL */

/* significant limbs of x, 0 if x is zero */
static mp_size_t
_norm(nn_srcptr x, mp_size_t n)
{
    while (n > 0 && x[n - 1] == 0)
        n--;
    return n;
}

static void
_mul_any_order(nn_ptr r, nn_srcptr a, mp_size_t an, nn_srcptr b, mp_size_t bn)
{
    if (an == 0 || bn == 0)
    {
        flint_mpn_zero(r, an + bn);
        return;
    }

    if (an >= bn)
        flint_mpn_mul(r, a, an, b, bn);
    else
        flint_mpn_mul(r, b, bn, a, an);
}

/*
    z receives (-1)^xs {x,xn} + (-1)^ys {y,yn} and the signed length is
    returned. The operand lengths are upper bounds, not normalized
    lengths: carrying an insignificant high limb through an addition
    costs a cycle, while normalizing every input costs a scan, and only
    the result needs to be canonical. z may be x or y exactly, and needs
    room for max(xn, yn) + 1 limbs.
*/
/* compare magnitudes given as upper-bound lengths, trimming both to
   their significant lengths; only the excess limbs are scanned, which is
   what makes it cheap enough to skip normalizing the inputs */
static int
_cmp_trim(nn_srcptr x, mp_size_t * xn, nn_srcptr y, mp_size_t * yn)
{
    while (*xn > *yn && x[*xn - 1] == 0)
        (*xn)--;
    while (*yn > *xn && y[*yn - 1] == 0)
        (*yn)--;

    if (*xn != *yn)
        return (*xn > *yn) ? 1 : -1;

    return (*xn == 0) ? 0 : mpn_cmp(x, y, *xn);
}

/* z == x permitted (y must be distinct). The manual's blanket rule for
   the mpn layer allows in-place operation "where source and destination
   are the same" and forbids partial overlap; for a source shorter than
   the destination but sharing its starting address -- which is what the
   swap branches of the general routine would pass as the second source
   of mpn_add/mpn_sub -- the text does not say which case applies. All
   implementations process low to high, so it is safe in practice; this
   variant expresses the aliasing through the equal-length primitives,
   whose in-place behavior is unambiguous, rather than relying on a
   reading of the rule. */
static slong
_signed_add_normalise_zx(nn_ptr z, mp_size_t xn, int xs,
                         nn_srcptr y, mp_size_t yn, int ys)
{
    mp_size_t rn;
    int rs;
    nn_srcptr x = z;

    /* the lengths are bounds: trim before treating an operand as
       absent, and never hand a zero length to the mpn primitives */
    while (xn > 0 && z[xn - 1] == 0)
        xn--;
    while (yn > 0 && y[yn - 1] == 0)
        yn--;
    if (xn == 0)
    {
        if (yn > 0)
            flint_mpn_copyi(z, y, yn);
        return ys ? -yn : yn;
    }
    if (yn == 0)
        return xs ? -xn : xn;

    if (xs == ys)
    {
        ulong cy;

        if (xn >= yn)
            cy = mpn_add(z, x, xn, y, yn);
        else
        {
            cy = mpn_add_n(z, x, y, xn);
            cy = mpn_add_1(z + xn, y + xn, yn - xn, cy);
        }
        rn = FLINT_MAX(xn, yn);
        if (cy != 0)
        {
            z[rn] = cy;
            rn++;
        }
        else
        {
            while (rn > 0 && z[rn - 1] == 0)
                rn--;
        }
        rs = xs;
    }
    else
    {
        int c = _cmp_trim(x, &xn, y, &yn);

        if (c == 0)
            return 0;

        if (c > 0)
        {
            mpn_sub(z, x, xn, y, yn);
            rn = xn;
            rs = xs;
        }
        else
        {
            ulong bw;

            bw = mpn_sub_n(z, y, x, xn);
            if (yn > xn)
                mpn_sub_1(z + xn, y + xn, yn - xn, bw);
            rn = yn;
            rs = ys;
        }

        while (rn > 0 && z[rn - 1] == 0)
            rn--;
        if (rn == 0)
            return 0;
    }

    return rs ? -rn : rn;
}

static slong
_signed_add_normalise(nn_ptr z, nn_srcptr x, mp_size_t xn, int xs,
                      nn_srcptr y, mp_size_t yn, int ys)
{
    mp_size_t rn;
    int rs;

    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    while (yn > 0 && y[yn - 1] == 0)
        yn--;
    if (xn == 0)
    {
        if (yn > 0)
            flint_mpn_copyi(z, y, yn);
        return ys ? -yn : yn;
    }
    if (yn == 0)
    {
        if (z != x)
            flint_mpn_copyi(z, x, xn);
        return xs ? -xn : xn;
    }

    if (xs == ys)
    {
        ulong cy;

        if (xn >= yn)
        {
            cy = mpn_add(z, x, xn, y, yn);
            rn = xn;
        }
        else
        {
            cy = mpn_add(z, y, yn, x, xn);
            rn = yn;
        }

        if (cy != 0)
        {
            z[rn] = cy;
            rn++;
        }
        else
        {
            /* the input lengths are bounds, not normalized: a product
               of trimmed operands may still carry a zero top limb, and
               the returned length must not include it (callers install
               it directly as an mpz size) */
            while (rn > 0 && z[rn - 1] == 0)
                rn--;
        }
        rs = xs;
    }
    else
    {
        int c = _cmp_trim(x, &xn, y, &yn);

        if (c == 0)
            return 0;

        if (c > 0)
        {
            mpn_sub(z, x, xn, y, yn);
            rn = xn;
            rs = xs;
        }
        else
        {
            mpn_sub(z, y, yn, x, xn);
            rn = yn;
            rs = ys;
        }
    }

    while (rn > 0 && z[rn - 1] == 0)
        rn--;

    return rs ? -rn : rn;
}

/*
    Is the Karatsuba scheme worth it for these shapes?

    It replaces the two cross products ar bi and ai br with the single
    auxiliary product (ar + ai)(br + bi), whose operands have
    max(arn, ain) + 1 and max(brn, bin) + 1 limbs. Comparing limb-product
    counts, that is a gain exactly when the auxiliary product is smaller
    than the two it removes -- true when the real and imaginary parts are
    internally balanced within each operand, false when any part is much
    shorter than its partner, without needing a ratio threshold. The two
    operands need not resemble each other.

    Verified directly on twelve unbalanced shapes (see the PR notes):
    eleven predictions correct, the miss on the exact boundary at a 4%
    cost -- a wrong call can only occur where the two algorithms are
    near-equal, since the counts being compared are what define the
    boundary.
*/
static int
_kara_worth_it(mp_size_t arn, mp_size_t ain, mp_size_t brn, mp_size_t bin,
               mp_size_t cutoff)
{
    mp_size_t maxa = FLINT_MAX(arn, ain);
    mp_size_t maxb = FLINT_MAX(brn, bin);

    if (maxa < cutoff || maxb < cutoff)
        return 0;

    /* the products below would overflow for parts of 2^31 limbs and
       more; every scheme is asymptotically better than schoolbook long
       before that */
    if ((maxa | maxb) >> (FLINT_BITS / 2 - 1) != 0)
        return 1;

    return (maxa + 1) * (maxb + 1) < arn * bin + ain * brn;
}

/*
    Schoolbook: four products, each output taking one of them directly
    and having the other combined into it in place. Nothing is normalized
    on the way in -- the products are formed at the lengths given and the
    additions run at those upper bounds.
*/
void
flint_mpn_mul_complex_classical(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC(FLINT_MAX(ain + bin, ain + brn) * sizeof(ulong));

    _mul_any_order(zr, ar, arn, br, brn);
    _mul_any_order(t, ai, ain, bi, bin);
    *zr_len = _signed_add_normalise_zx(zr, arn + brn, ar_sgn ^ br_sgn,
                                       t, ain + bin, ai_sgn ^ bi_sgn ^ 1);

    _mul_any_order(zi, ar, arn, bi, bin);
    _mul_any_order(t, ai, ain, br, brn);
    *zi_len = _signed_add_normalise_zx(zi, arn + bin, ar_sgn ^ bi_sgn,
                                       t, ain + brn, ai_sgn ^ br_sgn);

    TMP_END;
}

/*
    Karatsuba: t1 = ar br, t2 = ai bi, t3 = (ar + ai)(br + bi);
    zr = t1 - t2 and zi = t3 - (t1 + t2).

    Here the inputs *are* normalized first: an insignificant high limb
    costs a cycle in an addition but a whole extra row in a
    multiplication, and this path performs three of them.
*/
void
flint_mpn_mul_complex_karatsuba(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t maxa = FLINT_MAX(arn, ain), maxb = FLINT_MAX(brn, bin);
    nn_ptr t2, t3, u, v;
    mp_size_t un0 = maxa + 1, vn0 = maxb + 1;
    slong ul, vl, il;
    int s1, s2, s3;
    TMP_INIT;

    arn = _norm(ar, arn);
    ain = _norm(ai, ain);
    brn = _norm(br, brn);
    bin = _norm(bi, bin);

    /*
        The layout matters at large sizes, where fresh pages are a
        dominant linear cost: ar br is written straight into the
        caller's zr and both outputs are combined in place with the
        aliasing-tolerant variant above, so the scratch holds one
        product (ai bi), the product of the two operand sums, and the
        sums themselves -- about half the previous footprint -- and
        nothing is copied at the end. All combination values fit the
        caller's maxa + maxb + 1 window: the products are below
        2^(64(maxa+maxb)), the sum product below 2^(64(maxa+maxb)+2),
        and every intermediate below 2^(64(maxa+maxb)+3).
    */
    TMP_START;
    t2 = TMP_ALLOC((ain + bin + un0 + vn0 + un0 + vn0) * sizeof(ulong));
    t3 = t2 + ain + bin;
    u = t3 + un0 + vn0;
    v = u + un0;

    s1 = ar_sgn ^ br_sgn;
    s2 = ai_sgn ^ bi_sgn;

    ul = _signed_add_normalise(u, ar, arn, ar_sgn, ai, ain, ai_sgn);
    vl = _signed_add_normalise(v, br, brn, br_sgn, bi, bin, bi_sgn);
    _mul_any_order(t3, u, FLINT_ABS(ul), v, FLINT_ABS(vl));
    s3 = (ul < 0) ^ (vl < 0);

    _mul_any_order(t2, ai, ain, bi, bin);
    _mul_any_order(zr, ar, arn, br, brn);       /* t1, in place */

    /* zi = (t3 - t1) - t2; the first combination reads t1 from zr and
       aliases nothing, the second and the zr combination run in place */
    il = _signed_add_normalise(zi, t3, FLINT_ABS(ul) + FLINT_ABS(vl), s3,
                               zr, arn + brn, s1 ^ 1);
    *zi_len = _signed_add_normalise_zx(zi, FLINT_ABS(il), il < 0,
                                       t2, ain + bin, s2 ^ 1);

    *zr_len = _signed_add_normalise_zx(zr, arn + brn, s1,
                                       t2, ain + bin, s2 ^ 1);

    TMP_END;
}

/*
    Products with a zero part (a real or purely imaginary operand, or
    both): the transformed engine transforms only the nonzero parts and
    forms only the surviving part products, so the operand with one
    nonzero part is transformed once for both outputs (three forward and
    two inverse transforms for a real times a complex operand, against
    four and two for two independent real products). The cutoffs below
    send such shapes to it earlier than the general ones.
*/

slong flint_mpn_mul_complex_sparse_fft_cutoff = 300;
slong flint_mpn_mulmid_complex_sparse_fft_cutoff = 150;

/* whether some part is zero (the others nonzero or not) */
static int
_has_zero_part(nn_srcptr ar, mp_size_t arn, nn_srcptr ai, mp_size_t ain,
               nn_srcptr br, mp_size_t brn, nn_srcptr bi, mp_size_t bin)
{
    return _norm(ar, arn) == 0 || _norm(ai, ain) == 0 ||
           _norm(br, brn) == 0 || _norm(bi, bin) == 0;
}

/* the transformed path; returns 0 (leaving the outputs untouched) if
   fft_small is not built in or the plan is inadmissible */
int
flint_mpn_mul_complex_fft_small(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    mp_size_t maxa = FLINT_MAX(arn, ain), maxb = FLINT_MAX(brn, bin);

    return _fft_small_mpn_mul_complex(zr, zr_len, zi, zi_len, 0,
                maxa + maxb + 1, ar, arn, ar_sgn, ai, ain, ai_sgn,
                br, brn, br_sgn, bi, bin, bi_sgn);
#else
    return 0;
#endif
}

int
flint_mpn_sqr_complex_fft_small(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    mp_size_t maxa = FLINT_MAX(arn, ain);

    return _fft_small_mpn_sqr_complex(zr, zr_len, zi, zi_len, 0,
                2 * maxa + 1, ar, arn, ar_sgn, ai, ain, ai_sgn);
#else
    return 0;
#endif
}

static void
_mul_complex_general(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t maxa = FLINT_MAX(arn, ain), maxb = FLINT_MAX(brn, bin);

    /* a real or imaginary operand: the transformed path shares its
       transform between both output parts */
    if (FLINT_MIN(maxa, maxb) >= flint_mpn_mul_complex_sparse_fft_cutoff &&
        _has_zero_part(ar, arn, ai, ain, br, brn, bi, bin) &&
        flint_mpn_mul_complex_fft_small(zr, zr_len, zi, zi_len,
                ar, arn, ar_sgn, ai, ain, ai_sgn,
                br, brn, br_sgn, bi, bin, bi_sgn))
        return;

    if (_kara_worth_it(arn, ain, brn, bin, MUL_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (FLINT_MIN(maxa, maxb) >= flint_mpn_mul_complex_fft_cutoff &&
            flint_mpn_mul_complex_fft_small(zr, zr_len, zi, zi_len,
                    ar, arn, ar_sgn, ai, ain, ai_sgn,
                    br, brn, br_sgn, bi, bin, bi_sgn))
            return;

        flint_mpn_mul_complex_karatsuba(zr, zr_len, zi, zi_len, ar, arn,
                ar_sgn, ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
    }
    else
        flint_mpn_mul_complex_classical(zr, zr_len, zi, zi_len, ar, arn,
                ar_sgn, ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
}

/* zi = 2 ar ai for both squaring algorithms: the product written into
   zi and doubled in place */
static slong
_sqr_zi_double(nn_ptr zi, nn_srcptr ar, mp_size_t arn, int ar_sgn,
               nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t rn = arn + ain;
    ulong cy;

    _mul_any_order(zi, ar, arn, ai, ain);
    cy = mpn_lshift(zi, zi, rn, 1);
    zi[rn] = cy;
    rn += (cy != 0);
    while (rn > 0 && zi[rn - 1] == 0)
        rn--;

    return (rn == 0 || !(ar_sgn ^ ai_sgn)) ? rn : -rn;
}

void
flint_mpn_sqr_complex_classical(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    nn_ptr t1;
    TMP_INIT;

    TMP_START;
    t1 = TMP_ALLOC(2 * ain * sizeof(ulong));

    flint_mpn_sqr(zr, ar, arn);
    flint_mpn_sqr(t1, ai, ain);
    *zr_len = _signed_add_normalise_zx(zr, 2 * arn, 0, t1, 2 * ain, 1);

    *zi_len = _sqr_zi_double(zi, ar, arn, ar_sgn, ai, ain, ai_sgn);

    TMP_END;
}

void
flint_mpn_sqr_complex_karatsuba(nn_ptr zr, slong * zr_len,
    nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    /* zr = (ar + ai)(ar - ai): one multiplication where the classical
       form uses two squarings, so the inputs are normalized first */
    mp_size_t maxa = FLINT_MAX(arn, ain), rn;
    mp_size_t an = _norm(ar, arn), bn = _norm(ai, ain);
    nn_ptr u, v;
    slong ul, vl;
    TMP_INIT;

    TMP_START;
    u = TMP_ALLOC(2 * (maxa + 1) * sizeof(ulong));
    v = u + maxa + 1;

    /* exactly one of the two combinations is a magnitude sum (at most
       maxa + 1 limbs) and the other a magnitude difference (at most
       maxa), so their product spans at most 2 maxa + 1 limbs and is
       written straight into the caller's window; the scratch holds only
       the sums, and nothing is copied */
    ul = _signed_add_normalise(u, ar, an, ar_sgn, ai, bn, ai_sgn);
    vl = _signed_add_normalise(v, ar, an, ar_sgn, ai, bn, ai_sgn ^ 1);
    _mul_any_order(zr, u, FLINT_ABS(ul), v, FLINT_ABS(vl));
    rn = _norm(zr, FLINT_ABS(ul) + FLINT_ABS(vl));
    *zr_len = (rn == 0 || !((ul < 0) ^ (vl < 0))) ? rn : -rn;

    *zi_len = _sqr_zi_double(zi, ar, arn, ar_sgn, ai, ain, ai_sgn);

    TMP_END;
}

static void
_sqr_complex_general(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t maxa = FLINT_MAX(arn, ain);

    if (_kara_worth_it(arn, ain, arn, ain, SQR_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (maxa >= flint_mpn_sqr_complex_fft_cutoff &&
            flint_mpn_sqr_complex_fft_small(zr, zr_len, zi, zi_len,
                    ar, arn, ar_sgn, ai, ain, ai_sgn))
            return;

        flint_mpn_sqr_complex_karatsuba(zr, zr_len, zi, zi_len,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
    }
    else
        flint_mpn_sqr_complex_classical(zr, zr_len, zi, zi_len,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
}

/* z (n + 1 limbs) = (-1)^xs |x| + (-1)^ys |y| for n-limb magnitudes;
   returns the sign of the result */
static int
_signed_add_n1(nn_ptr z, nn_srcptr x, int xs, nn_srcptr y, int ys,
               mp_size_t n)
{
    if (xs == ys)
    {
        z[n] = mpn_add_n(z, x, y, n);
        return xs;
    }

    if (mpn_cmp(x, y, n) >= 0)
    {
        mpn_sub_n(z, x, y, n);
        z[n] = 0;
        return xs;
    }

    mpn_sub_n(z, y, x, n);
    z[n] = 0;
    return ys;
}

/*
    High variants: every product shares the length n, so the combinations
    are fixed-length signed adds -- one mpn_add_n, or one mpn_cmp plus
    mpn_sub_n -- in the style of nfloat's complex multiplication, with the
    carry limb completing the promised n + 1 and no accumulator to clear
    or copy out of.
*/
static void
_high_combine(nn_ptr z, int * zsgn, nn_srcptr x, int xs,
              nn_srcptr y, int ys, mp_size_t n)
{
    if (xs == ys)
    {
        z[n] = mpn_add_n(z, x, y, n);
        *zsgn = xs;
    }
    else if (mpn_cmp(x, y, n) >= 0)
    {
        mpn_sub_n(z, x, y, n);
        z[n] = 0;
        *zsgn = xs;
    }
    else
    {
        mpn_sub_n(z, y, x, n);
        z[n] = 0;
        *zsgn = ys;
    }
}

static void
_mul_complex_high(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn,
    nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn,
    nn_srcptr br, int br_sgn, nn_srcptr bi, int bi_sgn, mp_size_t n)
{
    TMP_INIT;
    TMP_START;

    if (n < MULHIGH_COMPLEX_KARATSUBA_CUTOFF)
    {
        nn_ptr h2 = TMP_ALLOC(2 * n * sizeof(ulong));
        nn_ptr h4 = h2 + n;

        flint_mpn_mulhigh_n(zr, ar, br, n);
        flint_mpn_mulhigh_n(h2, ai, bi, n);
        flint_mpn_mulhigh_n(zi, ar, bi, n);
        flint_mpn_mulhigh_n(h4, ai, br, n);

        _high_combine(zr, zr_sgn, zr, ar_sgn ^ br_sgn,
                      h2, ai_sgn ^ bi_sgn ^ 1, n);
        _high_combine(zi, zi_sgn, zi, ar_sgn ^ bi_sgn,
                      h4, ai_sgn ^ br_sgn, n);
    }
    else
    {
        /* three high products, the auxiliary one padded with a low zero
           limb so its window lines up with the others */
        /* zt holds h1 + h2 (n + 1 limbs) followed by the h3 - that
           combination (n + 3 limbs) */
        nn_ptr h1 = TMP_ALLOC((2 * n + (n + 2) + 2 * (n + 2)
                               + (2 * n + 4)) * sizeof(ulong));
        nn_ptr h2 = h1 + n, h3 = h2 + n;
        nn_ptr U = h3 + n + 2, V = U + n + 2, zt = V + n + 2;
        slong pl, il;
        int su, sv;

        flint_mpn_mulhigh_n(h1, ar, br, n);
        flint_mpn_mulhigh_n(h2, ai, bi, n);

        U[0] = 0;
        V[0] = 0;
        su = _signed_add_n1(U + 1, ar, ar_sgn, ai, ai_sgn, n);
        sv = _signed_add_n1(V + 1, br, br_sgn, bi, bi_sgn, n);
        flint_mpn_mulhigh_n(h3, U, V, n + 2);

        _high_combine(zr, zr_sgn, h1, ar_sgn ^ br_sgn,
                      h2, ai_sgn ^ bi_sgn ^ 1, n);

        /* zi = h3 - (h1 + h2), the sum in scratch since h3 spans n + 2 */
        pl = _signed_add_normalise(zt, h1, n, ar_sgn ^ br_sgn,
                                   h2, n, ai_sgn ^ bi_sgn);
        il = _signed_add_normalise(zt + n + 1, h3, n + 2, su ^ sv,
                                   zt, FLINT_ABS(pl), (pl < 0) ^ 1);
        {
            mp_size_t k = FLINT_MIN(FLINT_ABS(il), n + 1);

            flint_mpn_copyi(zi, zt + n + 1, k);
            if (k < n + 1)
                flint_mpn_zero(zi + k, n + 1 - k);
            *zi_sgn = (il < 0);
        }
    }

    TMP_END;
}

static void
_sqr_complex_high(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn,
    nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, mp_size_t n)
{
    TMP_INIT;
    TMP_START;

    if (n < SQRHIGH_COMPLEX_KARATSUBA_CUTOFF)
    {
        nn_ptr h1 = TMP_ALLOC(3 * n * sizeof(ulong));
        nn_ptr h2 = h1 + n, h3 = h2 + n;

        flint_mpn_sqrhigh_n(h1, ar, n);
        flint_mpn_sqrhigh_n(h2, ai, n);
        flint_mpn_mulhigh_n(h3, ar, ai, n);

        _high_combine(zr, zr_sgn, h1, 0, h2, 1, n);
        zi[n] = mpn_lshift(zi, h3, n, 1);
        *zi_sgn = ar_sgn ^ ai_sgn;
    }
    else
    {
        nn_ptr h1 = TMP_ALLOC(((n + 2) + n + 2 * (n + 2)) * sizeof(ulong));
        nn_ptr h2 = h1 + n + 2;
        nn_ptr U = h2 + n, V = U + n + 2;
        int su, sv;

        U[0] = 0;
        V[0] = 0;
        su = _signed_add_n1(U + 1, ar, ar_sgn, ai, ai_sgn, n);
        sv = _signed_add_n1(V + 1, ar, ar_sgn, ai, ai_sgn ^ 1, n);
        flint_mpn_mulhigh_n(h1, U, V, n + 2);
        flint_mpn_mulhigh_n(zi, ar, ai, n);

        flint_mpn_copyi(zr, h1, n + 1);
        *zr_sgn = su ^ sv;
        zi[n] = mpn_lshift(zi, zi, n, 1);
        *zi_sgn = ar_sgn ^ ai_sgn;
    }

    TMP_END;
}

void
flint_mpn_mul_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    FLINT_ASSERT(arn >= 1 && ain >= 1 && brn >= 1 && bin >= 1);

    _mul_complex_general(zr, zr_len, zi, zi_len, ar, arn, ar_sgn,
                         ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
}

void
flint_mpn_sqr_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    FLINT_ASSERT(arn >= 1 && ain >= 1);

    _sqr_complex_general(zr, zr_len, zi, zi_len, ar, arn, ar_sgn,
                         ai, ain, ai_sgn);
}

void
flint_mpn_mulhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn,
    nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn,
    nn_srcptr br, int br_sgn, nn_srcptr bi, int bi_sgn, mp_size_t n)
{
#if FLINT_HAVE_FFT_SMALL
    {
        slong lr, li;

        if (n >= flint_mpn_mul_complex_fft_cutoff &&
            _fft_small_mpn_mul_complex(zr, &lr, zi, &li, n, 2 * n + 1,
                ar, n, ar_sgn, ai, n, ai_sgn, br, n, br_sgn, bi, n, bi_sgn))
        {
            _high_pad(zr, zr_sgn, lr, n + 1);
            _high_pad(zi, zi_sgn, li, n + 1);
            return;
        }
    }
#endif
    _mul_complex_high(zr, zr_sgn, zi, zi_sgn, ar, ar_sgn, ai, ai_sgn,
                      br, br_sgn, bi, bi_sgn, n);
}

void
flint_mpn_sqrhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn,
    nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, mp_size_t n)
{
#if FLINT_HAVE_FFT_SMALL
    {
        slong lr, li;

        if (n >= flint_mpn_sqr_complex_fft_cutoff &&
            _fft_small_mpn_sqr_complex(zr, &lr, zi, &li, n, 2 * n + 1,
                ar, n, ar_sgn, ai, n, ai_sgn))
        {
            _high_pad(zr, zr_sgn, lr, n + 1);
            _high_pad(zi, zi_sgn, li, n + 1);
            return;
        }
    }
#endif
    _sqr_complex_high(zr, zr_sgn, zi, zi_sgn, ar, ar_sgn, ai, ai_sgn, n);
}

/*
    Low products: zr and zi receive the real and imaginary parts of the
    product mod B^n as n-limb two's-complement residues. The inputs are
    magnitudes with sign bits as for the full products, of any lengths
    (>= 1 limb, possibly beyond n); only their low n limbs matter, since
    (-1)^s |x| = (-1)^s (|x| mod B^n) mod B^n. The residues need no sign:
    they are the low limbs of the exact parts in two's complement, which
    is what a Kronecker unpacking of signed fields reads.
*/

slong flint_mpn_mullow_complex_fft_cutoff = 230;
slong flint_mpn_sqrlow_complex_fft_cutoff = 270;

/* t (n limbs) = |a| |b| mod B^n; zero lengths allowed */
static void
_mullow_any(nn_ptr t, mp_size_t n, nn_srcptr a, mp_size_t an,
            nn_srcptr b, mp_size_t bn)
{
    mp_size_t zh;

    an = _norm(a, FLINT_MIN(an, n));
    bn = _norm(b, FLINT_MIN(bn, n));

    if (an == 0 || bn == 0)
    {
        flint_mpn_zero(t, n);
        return;
    }

    zh = FLINT_MIN(n, an + bn);

    if (an >= bn)
        flint_mpn_mulmid(t, a, an, b, bn, 0, zh);
    else
        flint_mpn_mulmid(t, b, bn, a, an, 0, zh);

    if (zh < n)
        flint_mpn_zero(t + zh, n - zh);
}

/* z = (-1)^xs x + (-1)^ys y mod B^n; z may alias x or y */
static void
_addsub_mod(nn_ptr z, nn_srcptr x, int xs, nn_srcptr y, int ys, mp_size_t n)
{
    if (xs == ys)
    {
        mpn_add_n(z, x, y, n);
        if (xs)
            mpn_neg(z, z, n);
    }
    else if (xs)
        mpn_sub_n(z, y, x, n);
    else
        mpn_sub_n(z, x, y, n);
}

/* z = 2 x mod B^n, negated if s */
static void
_double_mod(nn_ptr z, nn_srcptr x, int s, mp_size_t n)
{
    mpn_lshift(z, x, n, 1);
    if (s)
        mpn_neg(z, z, n);
}

void
flint_mpn_mullow_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC(n * sizeof(ulong));

    _mullow_any(zr, n, ar, arn, br, brn);
    _mullow_any(t, n, ai, ain, bi, bin);
    _addsub_mod(zr, zr, ar_sgn ^ br_sgn, t, ai_sgn ^ bi_sgn ^ 1, n);

    _mullow_any(zi, n, ar, arn, bi, bin);
    _mullow_any(t, n, ai, ain, br, brn);
    _addsub_mod(zi, zi, ar_sgn ^ bi_sgn, t, ai_sgn ^ br_sgn, n);

    TMP_END;
}

void
flint_mpn_mullow_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t maxa, maxb;
    nn_ptr t2, u, v;
    slong ul, vl;
    int s1, s2, s3;
    TMP_INIT;

    arn = _norm(ar, FLINT_MIN(arn, n));
    ain = _norm(ai, FLINT_MIN(ain, n));
    brn = _norm(br, FLINT_MIN(brn, n));
    bin = _norm(bi, FLINT_MIN(bin, n));
    maxa = FLINT_MAX(arn, ain);
    maxb = FLINT_MAX(brn, bin);

    TMP_START;
    t2 = TMP_ALLOC((n + maxa + 1 + maxb + 1) * sizeof(ulong));
    u = t2 + n;
    v = u + maxa + 1;

    s1 = ar_sgn ^ br_sgn;
    s2 = ai_sgn ^ bi_sgn;

    /* the sums of the truncated parts are congruent to the sums of the
       parts, and so is their truncated product */
    ul = _signed_add_normalise(u, ar, arn, ar_sgn, ai, ain, ai_sgn);
    vl = _signed_add_normalise(v, br, brn, br_sgn, bi, bin, bi_sgn);
    s3 = (ul < 0) ^ (vl < 0);

    _mullow_any(zi, n, u, FLINT_ABS(ul), v, FLINT_ABS(vl));    /* t3 */
    _mullow_any(t2, n, ai, ain, bi, bin);
    _mullow_any(zr, n, ar, arn, br, brn);                      /* t1 */

    /* zi = t3 - t1 - t2, zr = t1 - t2 */
    _addsub_mod(zi, zi, s3, zr, s1 ^ 1, n);
    _addsub_mod(zi, zi, 0, t2, s2 ^ 1, n);
    _addsub_mod(zr, zr, s1, t2, s2 ^ 1, n);

    TMP_END;
}

void
flint_mpn_sqrlow_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    nn_ptr t;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC(n * sizeof(ulong));

    _mullow_any(zr, n, ar, arn, ar, arn);
    _mullow_any(t, n, ai, ain, ai, ain);
    mpn_sub_n(zr, zr, t, n);

    _mullow_any(t, n, ar, arn, ai, ain);
    _double_mod(zi, t, ar_sgn ^ ai_sgn, n);

    TMP_END;
}

void
flint_mpn_sqrlow_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t maxa;
    nn_ptr u, v;
    slong ul, vl;
    TMP_INIT;

    arn = _norm(ar, FLINT_MIN(arn, n));
    ain = _norm(ai, FLINT_MIN(ain, n));
    maxa = FLINT_MAX(arn, ain);

    TMP_START;
    u = TMP_ALLOC(2 * (maxa + 1) * sizeof(ulong));
    v = u + maxa + 1;

    /* zr = (ar + ai)(ar - ai) */
    ul = _signed_add_normalise(u, ar, arn, ar_sgn, ai, ain, ai_sgn);
    vl = _signed_add_normalise(v, ar, arn, ar_sgn, ai, ain, ai_sgn ^ 1);
    _mullow_any(zr, n, u, FLINT_ABS(ul), v, FLINT_ABS(vl));
    if ((ul < 0) ^ (vl < 0))
        mpn_neg(zr, zr, n);

    _mullow_any(zi, n, ar, arn, ai, ain);
    _double_mod(zi, zi, ar_sgn ^ ai_sgn, n);

    TMP_END;
}

int
flint_mpn_mullow_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    /* the residues mod B^n are the windows [0, n) */
    return _fft_small_mpn_mulmid_complex(zr, zi, 0, n, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
#else
    return 0;
#endif
}

int
flint_mpn_sqrlow_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    return _fft_small_mpn_sqrmid_complex(zr, zi, 0, n, ar, arn, ar_sgn,
                                         ai, ain, ai_sgn);
#else
    return 0;
#endif
}

void
flint_mpn_mullow_complex(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t maxa, maxb;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(arn >= 1 && ain >= 1 && brn >= 1 && bin >= 1);

    /* only the low n limbs of each part matter */
    arn = FLINT_MIN(arn, n);
    ain = FLINT_MIN(ain, n);
    brn = FLINT_MIN(brn, n);
    bin = FLINT_MIN(bin, n);
    maxa = FLINT_MAX(arn, ain);
    maxb = FLINT_MAX(brn, bin);

    if (FLINT_MIN(maxa, maxb) >= flint_mpn_mul_complex_sparse_fft_cutoff &&
        _has_zero_part(ar, arn, ai, ain, br, brn, bi, bin) &&
        flint_mpn_mullow_complex_fft_small(zr, zi, n,
                ar, arn, ar_sgn, ai, ain, ai_sgn,
                br, brn, br_sgn, bi, bin, bi_sgn))
        return;

    if (_kara_worth_it(arn, ain, brn, bin, MUL_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (FLINT_MIN(maxa, maxb) >= flint_mpn_mullow_complex_fft_cutoff &&
            flint_mpn_mullow_complex_fft_small(zr, zi, n,
                    ar, arn, ar_sgn, ai, ain, ai_sgn,
                    br, brn, br_sgn, bi, bin, bi_sgn))
            return;

        flint_mpn_mullow_complex_karatsuba(zr, zi, n, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
    }
    else
        flint_mpn_mullow_complex_classical(zr, zi, n, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
}

void
flint_mpn_sqrlow_complex(nn_ptr zr, nn_ptr zi, mp_size_t n,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t maxa;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(arn >= 1 && ain >= 1);

    arn = FLINT_MIN(arn, n);
    ain = FLINT_MIN(ain, n);
    maxa = FLINT_MAX(arn, ain);

    if (_kara_worth_it(arn, ain, arn, ain, SQR_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (maxa >= flint_mpn_sqrlow_complex_fft_cutoff &&
            flint_mpn_sqrlow_complex_fft_small(zr, zi, n,
                    ar, arn, ar_sgn, ai, ain, ai_sgn))
            return;

        flint_mpn_sqrlow_complex_karatsuba(zr, zi, n,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
    }
    else
        flint_mpn_sqrlow_complex_classical(zr, zi, n,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
}

/*
    Windowed products: zr and zi receive the limbs [lo, hi) of the real
    and imaginary parts in two's complement, i.e. floor(x / B^lo) + e
    mod B^(hi - lo) for each exact part x, with a small error e in the
    lowest limb (none for lo = 0, where these are the low products).
    The schoolbook and Karatsuba variants combine real windows from
    flint_mpn_mulmid, each a lower approximation short by less than
    min(an, bn) B units; the transformed variant has |e| <= 1.
*/

slong flint_mpn_mulmid_complex_fft_cutoff = 120;
slong flint_mpn_sqrmid_complex_fft_cutoff = 160;

/* t (hi - lo limbs) = limbs [lo, hi) of |a| |b|, as flint_mpn_mulmid
   computes them (a lower approximation); any lengths */
static void
_mulmid_any(nn_ptr t, mp_size_t lo, mp_size_t hi, nn_srcptr a, mp_size_t an,
            nn_srcptr b, mp_size_t bn)
{
    mp_size_t zh;

    an = _norm(a, FLINT_MIN(an, hi));
    bn = _norm(b, FLINT_MIN(bn, hi));

    if (an == 0 || bn == 0 || lo >= an + bn)
    {
        flint_mpn_zero(t, hi - lo);
        return;
    }

    zh = FLINT_MIN(hi, an + bn);

    if (an >= bn)
        flint_mpn_mulmid(t, a, an, b, bn, lo, zh);
    else
        flint_mpn_mulmid(t, b, bn, a, an, lo, zh);

    if (zh < hi)
        flint_mpn_zero(t + zh - lo, hi - zh);
}

void
flint_mpn_mulmid_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t w = hi - lo;
    nn_ptr t;
    TMP_INIT;

    if (lo == 0)
    {
        flint_mpn_mullow_complex_classical(zr, zi, hi, ar, arn, ar_sgn,
            ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
        return;
    }

    TMP_START;
    t = TMP_ALLOC(w * sizeof(ulong));

    _mulmid_any(zr, lo, hi, ar, arn, br, brn);
    _mulmid_any(t, lo, hi, ai, ain, bi, bin);
    _addsub_mod(zr, zr, ar_sgn ^ br_sgn, t, ai_sgn ^ bi_sgn ^ 1, w);

    _mulmid_any(zi, lo, hi, ar, arn, bi, bin);
    _mulmid_any(t, lo, hi, ai, ain, br, brn);
    _addsub_mod(zi, zi, ar_sgn ^ bi_sgn, t, ai_sgn ^ br_sgn, w);

    TMP_END;
}

void
flint_mpn_mulmid_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t w = hi - lo, maxa, maxb;
    nn_ptr t2, u, v;
    slong ul, vl;
    int s1, s2, s3;
    TMP_INIT;

    if (lo == 0)
    {
        flint_mpn_mullow_complex_karatsuba(zr, zi, hi, ar, arn, ar_sgn,
            ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
        return;
    }

    arn = _norm(ar, FLINT_MIN(arn, hi));
    ain = _norm(ai, FLINT_MIN(ain, hi));
    brn = _norm(br, FLINT_MIN(brn, hi));
    bin = _norm(bi, FLINT_MIN(bin, hi));
    maxa = FLINT_MAX(arn, ain);
    maxb = FLINT_MAX(brn, bin);

    TMP_START;
    t2 = TMP_ALLOC((w + maxa + 1 + maxb + 1) * sizeof(ulong));
    u = t2 + w;
    v = u + maxa + 1;

    s1 = ar_sgn ^ br_sgn;
    s2 = ai_sgn ^ bi_sgn;

    ul = _signed_add_normalise(u, ar, arn, ar_sgn, ai, ain, ai_sgn);
    vl = _signed_add_normalise(v, br, brn, br_sgn, bi, bin, bi_sgn);
    s3 = (ul < 0) ^ (vl < 0);

    _mulmid_any(zi, lo, hi, u, FLINT_ABS(ul), v, FLINT_ABS(vl));    /* t3 */
    _mulmid_any(t2, lo, hi, ai, ain, bi, bin);
    _mulmid_any(zr, lo, hi, ar, arn, br, brn);                      /* t1 */

    _addsub_mod(zi, zi, s3, zr, s1 ^ 1, w);
    _addsub_mod(zi, zi, 0, t2, s2 ^ 1, w);
    _addsub_mod(zr, zr, s1, t2, s2 ^ 1, w);

    TMP_END;
}

void
flint_mpn_sqrmid_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t w = hi - lo;
    nn_ptr t;
    TMP_INIT;

    if (lo == 0)
    {
        flint_mpn_sqrlow_complex_classical(zr, zi, hi, ar, arn, ar_sgn, ai, ain, ai_sgn);
        return;
    }

    TMP_START;
    t = TMP_ALLOC(w * sizeof(ulong));

    _mulmid_any(zr, lo, hi, ar, arn, ar, arn);
    _mulmid_any(t, lo, hi, ai, ain, ai, ain);
    mpn_sub_n(zr, zr, t, w);

    /* 2 ar ai from the window one limb lower, so that the doubling
       shifts in the true bit below lo rather than a zero */
    {
        nn_ptr t2 = TMP_ALLOC((w + 1) * sizeof(ulong));
        _mulmid_any(t2, lo - 1, hi, ar, arn, ai, ain);
        mpn_lshift(t2, t2, w + 1, 1);
        flint_mpn_copyi(zi, t2 + 1, w);
        if (ar_sgn ^ ai_sgn)
            mpn_neg(zi, zi, w);
    }

    TMP_END;
}

void
flint_mpn_sqrmid_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t w = hi - lo, maxa;
    nn_ptr u, v, t2;
    slong ul, vl;
    TMP_INIT;

    if (lo == 0)
    {
        flint_mpn_sqrlow_complex_karatsuba(zr, zi, hi, ar, arn, ar_sgn, ai, ain, ai_sgn);
        return;
    }

    arn = _norm(ar, FLINT_MIN(arn, hi));
    ain = _norm(ai, FLINT_MIN(ain, hi));
    maxa = FLINT_MAX(arn, ain);

    TMP_START;
    u = TMP_ALLOC((2 * (maxa + 1) + w + 1) * sizeof(ulong));
    v = u + maxa + 1;
    t2 = v + maxa + 1;

    /* zr = (ar + ai)(ar - ai) */
    ul = _signed_add_normalise(u, ar, arn, ar_sgn, ai, ain, ai_sgn);
    vl = _signed_add_normalise(v, ar, arn, ar_sgn, ai, ain, ai_sgn ^ 1);
    _mulmid_any(zr, lo, hi, u, FLINT_ABS(ul), v, FLINT_ABS(vl));
    if ((ul < 0) ^ (vl < 0))
        mpn_neg(zr, zr, w);

    _mulmid_any(t2, lo - 1, hi, ar, arn, ai, ain);
    mpn_lshift(t2, t2, w + 1, 1);
    flint_mpn_copyi(zi, t2 + 1, w);
    if (ar_sgn ^ ai_sgn)
        mpn_neg(zi, zi, w);

    TMP_END;
}

int
flint_mpn_mulmid_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    return _fft_small_mpn_mulmid_complex(zr, zi, lo, hi, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
#else
    return 0;
#endif
}

int
flint_mpn_sqrmid_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
#if FLINT_HAVE_FFT_SMALL
    return _fft_small_mpn_sqrmid_complex(zr, zi, lo, hi, ar, arn, ar_sgn,
                                         ai, ain, ai_sgn);
#else
    return 0;
#endif
}

void
flint_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn)
{
    mp_size_t maxa, maxb;

    FLINT_ASSERT(0 <= lo && lo < hi);
    FLINT_ASSERT(arn >= 1 && ain >= 1 && brn >= 1 && bin >= 1);

    if (lo == 0)
    {
        flint_mpn_mullow_complex(zr, zi, hi, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                 br, brn, br_sgn, bi, bin, bi_sgn);
        return;
    }

    arn = FLINT_MIN(arn, hi);
    ain = FLINT_MIN(ain, hi);
    brn = FLINT_MIN(brn, hi);
    bin = FLINT_MIN(bin, hi);
    maxa = FLINT_MAX(arn, ain);
    maxb = FLINT_MAX(brn, bin);

    if (FLINT_MIN(maxa, maxb) >= flint_mpn_mulmid_complex_sparse_fft_cutoff &&
        _has_zero_part(ar, arn, ai, ain, br, brn, bi, bin) &&
        flint_mpn_mulmid_complex_fft_small(zr, zi, lo, hi,
                ar, arn, ar_sgn, ai, ain, ai_sgn,
                br, brn, br_sgn, bi, bin, bi_sgn))
        return;

    if (_kara_worth_it(arn, ain, brn, bin, MUL_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (FLINT_MIN(maxa, maxb) >= flint_mpn_mulmid_complex_fft_cutoff &&
            flint_mpn_mulmid_complex_fft_small(zr, zi, lo, hi,
                    ar, arn, ar_sgn, ai, ain, ai_sgn,
                    br, brn, br_sgn, bi, bin, bi_sgn))
            return;

        flint_mpn_mulmid_complex_karatsuba(zr, zi, lo, hi, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
    }
    else
        flint_mpn_mulmid_complex_classical(zr, zi, lo, hi, ar, arn, ar_sgn,
                ai, ain, ai_sgn, br, brn, br_sgn, bi, bin, bi_sgn);
}

void
flint_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn)
{
    mp_size_t maxa;

    FLINT_ASSERT(0 <= lo && lo < hi);
    FLINT_ASSERT(arn >= 1 && ain >= 1);

    if (lo == 0)
    {
        flint_mpn_sqrlow_complex(zr, zi, hi, ar, arn, ar_sgn, ai, ain, ai_sgn);
        return;
    }

    arn = FLINT_MIN(arn, hi);
    ain = FLINT_MIN(ain, hi);
    maxa = FLINT_MAX(arn, ain);

    if (_kara_worth_it(arn, ain, arn, ain, SQR_COMPLEX_KARATSUBA_CUTOFF))
    {
        if (maxa >= flint_mpn_sqrmid_complex_fft_cutoff &&
            flint_mpn_sqrmid_complex_fft_small(zr, zi, lo, hi,
                    ar, arn, ar_sgn, ai, ain, ai_sgn))
            return;

        flint_mpn_sqrmid_complex_karatsuba(zr, zi, lo, hi,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
    }
    else
        flint_mpn_sqrmid_complex_classical(zr, zi, lo, hi,
                ar, arn, ar_sgn, ai, ain, ai_sgn);
}

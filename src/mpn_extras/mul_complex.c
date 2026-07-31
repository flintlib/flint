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
#include "gr.h"
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
    need no case analysis. Above the threshold the operands are converted
    once into a transformed-integer context and reused across the
    pointwise products: 4 forward + 2 inverse transforms replace the
    8 + 4 of four fused multiplications (2 + 2 replace 6 + 3 for
    squaring).

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

slong flint_mpn_mul_complex_fft_cutoff = 700;
slong flint_mpn_sqr_complex_fft_cutoff = 3000000;

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

/* convert an accumulated element out into (z, zlimbs) starting at limb
   offset 'shift' of the magnitude; returns 0 on any refusal. shift > 0
   uses the truncated conversion (at most one unit of error in the
   lowest returned limb, the mulhigh contract); shift = 0 is exact. */
static int
_tmpn_out(nn_ptr z, mp_size_t zlimbs, slong * zlen, gr_ptr acc,
          mp_size_t shift, nn_ptr t, ulong tmax, gr_ctx_t tctx)
{
    slong need, tn;
    int sg, status;

    /* the accumulators are dead after this, so the conversion may
       consume them and skip copying the transform; t is a slice of a
       dead operand's storage, large enough by construction (an element
       holds np * stride doubles against need ~ zlimbs + a few limbs) */
    if (shift > 0)
    {
        need = gr_transformed_mpn_get_limbs_trunc(tctx, acc, shift);
        if ((ulong) need > tmax)
            return 0;
        status = gr_transformed_mpn_get_trunc_destructive(t, need, &tn, &sg,
                                                          shift, acc, tctx);
    }
    else
    {
        need = gr_transformed_mpn_get_limbs(tctx, acc);
        if ((ulong) need > tmax)
            return 0;
        status = gr_transformed_mpn_get_destructive(t, need, &tn, &sg, acc,
                                                    tctx);
    }

    if (status == GR_SUCCESS)
    {
        /* the true value must fit the caller's window */
        if (tn > zlimbs)
            status = GR_UNABLE;
        else
        {
            if (tn > 0)
                flint_mpn_copyi(z, t, tn);
            *zlen = (tn == 0 || !sg) ? tn : -(slong) tn;
        }
    }
    return status == GR_SUCCESS;
}

/* fft path shared by mul and mulhigh: shift = 0 writes 2n+1 limbs,
   shift = n writes n+1 limbs */
static int
_mul_complex_fft(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    nn_srcptr br, mp_size_t brn, int br_sgn,
    nn_srcptr bi, mp_size_t bin, int bi_sgn,
    mp_size_t shift, mp_size_t zlimbs)
{
    gr_ctx_t tctx;
    gr_ptr E;
    double * base;
    int ok = 1;
#define E_(i) GR_ENTRY(E, i, tctx->sizeof_elem)

    if (gr_ctx_init_transformed_mpn(tctx,
            FLINT_BITS * (slong) (FLINT_MAX(arn, ain) + FLINT_MAX(brn, bin))
            + 8, 2, 1, 6) != GR_SUCCESS)
        return 0;

    {
        ulong esz = gr_transformed_mpn_sizeof_data(tctx);

        base = (double *) mpn_ctx_fit_buffer(get_default_mpn_ctx(),
                    6 * esz + 6 * tctx->sizeof_elem);
        E = (gr_ptr) (base + 6 * (esz / sizeof(double)));
        { slong i; for (i = 0; i < 6; i++)
            gr_transformed_mpn_init_borrowed(E_(i),
                    base + i * (esz / sizeof(double)), tctx); }
    }

    ok = ok && gr_transformed_mpn_set(E_(0), ar, arn, ar_sgn, tctx) == GR_SUCCESS;
    ok = ok && gr_transformed_mpn_set(E_(1), ai, ain, ai_sgn, tctx) == GR_SUCCESS;
    ok = ok && gr_transformed_mpn_set(E_(2), br, brn, br_sgn, tctx) == GR_SUCCESS;
    ok = ok && gr_transformed_mpn_set(E_(3), bi, bin, bi_sgn, tctx) == GR_SUCCESS;

    /* zr = ar br - ai bi, zi = ar bi + ai br */
    ok = ok && gr_mul(E_(4), E_(0), E_(2), tctx) == GR_SUCCESS;
    ok = ok && gr_submul(E_(4), E_(1), E_(3), tctx) == GR_SUCCESS;
    ok = ok && gr_mul(E_(5), E_(0), E_(3), tctx) == GR_SUCCESS;
    ok = ok && gr_addmul(E_(5), E_(1), E_(2), tctx) == GR_SUCCESS;

    /* the operands are dead once the pointwise stage is done: their
       storage (the first two scratch slices) doubles as the export
       staging buffers */
    {
        ulong esz = gr_transformed_mpn_sizeof_data(tctx);
        ulong tmax = esz / sizeof(ulong);

        ok = ok && _tmpn_out(zr, zlimbs, zr_len, E_(4), shift,
                             (nn_ptr) base, tmax, tctx);
        ok = ok && _tmpn_out(zi, zlimbs, zi_len, E_(5), shift,
                             (nn_ptr) (base + esz / sizeof(double)), tmax,
                             tctx);
    }

    { slong i; for (i = 0; i < 6; i++) gr_clear(E_(i), tctx); }
    gr_ctx_clear(tctx);
#undef E_
    return ok;
}

static int
_sqr_complex_fft(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len,
    nn_srcptr ar, mp_size_t arn, int ar_sgn,
    nn_srcptr ai, mp_size_t ain, int ai_sgn,
    mp_size_t shift, mp_size_t zlimbs)
{
    gr_ctx_t tctx;
    gr_ptr E;
    double * base;
    int ok = 1;
#define E_(i) GR_ENTRY(E, i, tctx->sizeof_elem)

    /*
        terms_bound is 2, not 4: an earlier version formed
        (ar + ai)(ar - ai) on sums of transforms, whose slot values need
        four products' capacity -- the profile then packs fewer bits per
        slot, and since the conversions dominate (their cost is
        proportional to the slot count), both directions inflated by more
        than the two forward transforms the method saves over the
        two-multiplication Karatsuba square. The identity buys nothing
        where multiplications are pointwise; the direct
        ar * ar - ai * ai keeps the same geometry as the complex product.
    */
    if (gr_ctx_init_transformed_mpn(tctx,
            FLINT_BITS * (slong) (2 * FLINT_MAX(arn, ain)) + 8, 2, 1, 4)
            != GR_SUCCESS)
        return 0;

    {
        ulong esz = gr_transformed_mpn_sizeof_data(tctx);

        base = (double *) mpn_ctx_fit_buffer(get_default_mpn_ctx(),
                    4 * esz + 4 * tctx->sizeof_elem);
        E = (gr_ptr) (base + 4 * (esz / sizeof(double)));
        { slong i; for (i = 0; i < 4; i++)
            gr_transformed_mpn_init_borrowed(E_(i),
                    base + i * (esz / sizeof(double)), tctx); }
    }

    /* the operands enter as magnitudes: the squares are sign free and
       2 ar ai has the known sign ar_sgn ^ ai_sgn, applied below -- so
       the imaginary part keeps a nonnegative element and takes the
       cheaper unsigned reconstruction, which lacks the bias and
       two's-complement sign-resolution passes of the signed one; only
       the real part, whose sign genuinely depends on the values, pays
       for the signed export */
    ok = ok && gr_transformed_mpn_set(E_(0), ar, arn, 0, tctx) == GR_SUCCESS;
    ok = ok && gr_transformed_mpn_set(E_(1), ai, ain, 0, tctx) == GR_SUCCESS;

    /* zr = ar ar - ai ai; zi = ar ai doubled by a pointwise addition,
       which costs a pass over the evaluations against a full pointwise
       multiplication */
    ok = ok && gr_mul(E_(2), E_(0), E_(0), tctx) == GR_SUCCESS;
    ok = ok && gr_submul(E_(2), E_(1), E_(1), tctx) == GR_SUCCESS;
    ok = ok && gr_mul(E_(3), E_(0), E_(1), tctx) == GR_SUCCESS;
    ok = ok && gr_add(E_(3), E_(3), E_(3), tctx) == GR_SUCCESS;

    /* operands dead after the pointwise stage; their slices stage the
       exports */
    {
        ulong esz = gr_transformed_mpn_sizeof_data(tctx);
        ulong tmax = esz / sizeof(ulong);

        ok = ok && _tmpn_out(zr, zlimbs, zr_len, E_(2), shift,
                             (nn_ptr) base, tmax, tctx);
        ok = ok && _tmpn_out(zi, zlimbs, zi_len, E_(3), shift,
                             (nn_ptr) (base + esz / sizeof(double)), tmax,
                             tctx);

        /* attach the known sign of 2 ar ai */
        if (ok && (ar_sgn ^ ai_sgn))
            *zi_len = -*zi_len;
    }

    { slong i; for (i = 0; i < 4; i++) gr_clear(E_(i), tctx); }
    gr_ctx_clear(tctx);
#undef E_
    return ok;
}

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

    return _mul_complex_fft(zr, zr_len, zi, zi_len, ar, arn, ar_sgn,
                            ai, ain, ai_sgn, br, brn, br_sgn, bi, bin,
                            bi_sgn, 0, maxa + maxb + 1);
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

    return _sqr_complex_fft(zr, zr_len, zi, zi_len, ar, arn, ar_sgn,
                            ai, ain, ai_sgn, 0, 2 * maxa + 1);
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

        flint_mpn_sqrhigh(h1, ar, n);
        flint_mpn_sqrhigh(h2, ai, n);
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
            _mul_complex_fft(zr, &lr, zi, &li, ar, n, ar_sgn, ai, n, ai_sgn,
                             br, n, br_sgn, bi, n, bi_sgn, n, n + 1))
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
            _sqr_complex_fft(zr, &lr, zi, &li, ar, n, ar_sgn, ai, n, ai_sgn,
                             n, n + 1))
        {
            _high_pad(zr, zr_sgn, lr, n + 1);
            _high_pad(zi, zi_sgn, li, n + 1);
            return;
        }
    }
#endif
    _sqr_complex_high(zr, zr_sgn, zi, zi_sgn, ar, ar_sgn, ai, ai_sgn, n);
}

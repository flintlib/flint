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
#include "fft_small.h"
#endif

/*
    General high products: the top n limbs of x * y (or x^2) for operands of
    any lengths, together with a guard limb (the limb just below), as a
    lower approximation with the deficit bounds FLINT_MPN_MULHIGH_GUARD_ERR_BOUND
    and friends.

    Writing zn = xn + yn, the product has zn limbs and the guard limb is
    limb d = zn - n - 1 of it. Operands longer than n limbs are first
    truncated to their top n limbs: dropping the low limbs of x loses
    x_lo * y < B^(xn-n) * B^yn = B^(d+1), that is, less than B units of the
    guard limb, of which the dominant part, the high half of the product of
    the limb just below the cut with the top limb of y, is added back. What
    then remains of the dropped part is less than 3 units of the guard limb
    per truncated operand (the low half of that product, the limb below the
    cut times the rest of y, and everything further down each contribute
    less than one unit).

    With both operands of at most n limbs the shapes are handled as
    follows (res = n limbs, zn = xn + yn, d = zn - n - 1):

      xn = yn = n           flint_mpn_mulhigh_n (deficit <= 2n - 4)
      zn <= n               the full product, zero padded above (exact)
      small n               the full product on the stack (cheaper than
                            zero padding the operands for the n x n kernel)
      large min(xn, yn)     the window [d - 1, zn) through fft_small, with
                            the two lowest limbs diverted to a separate
                            buffer (the lowest absorbs the carry deficit of
                            the window, the next is the guard limb)
      min(xn, yn) >= 7n/8   operands zero padded below to n limbs for the
                            balanced high product (same limbs of the
                            product, same deficit)
      xn = n > yn           x = x_hi B^yn + x_lo: the full product x_hi * y
                            lands exactly on the output, plus the balanced
                            high product of x_lo * y (deficit <= 2yn - 4)
      xn, yn < n            the full product (splitting off high products
                            does not pay in the Karatsuba range)

    so that the deficit before truncation never exceeds 2n - 4 (the window
    deficit of fft_small is below min(xn, yn) units of the guard limb) and
    the total is below (2n - 3) + 6, i.e. at most 2n + 2 as an integer for
    n >= 3, and at most 5 for n <= 2 where the kernels are exact. The
    macros state 2n + 4.
*/

/* a * b with the operands in either order */
static void
_mul_any(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    if (an >= bn)
        flint_mpn_mul(r, a, an, b, bn);
    else
        flint_mpn_mul(r, b, bn, a, an);
}

/* Below this many output limbs a full product (with its own kernel)
   followed by a copy beats every split: the balanced high product kernels
   only save their fraction of the work on larger operands. Measured on
   Zen 3 (where padding or peeling wins from n = 8) and Skylake. */
#define MULHIGH_FULL_PRODUCT_MAX_N 6

/* Max product length computed on the stack in the small range. */
#define MULHIGH_STACK_LIMBS (2 * FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH + 2)

/* Operands of at least this fraction of n limbs are zero padded to n limbs
   for the balanced high product; below it, an operand of full length is
   split (x_hi * y in full plus the high product of x_lo * y) and two
   short operands are multiplied in full. With the balanced high product
   at about 0.65 of the cost of the full product of the same size and
   Karatsuba scaling, padding both operands from m to n limbs costs
   0.65 (n/m)^1.58 relative to the full m x m product, a win from about
   m >= 0.8 n; measured (Zen 3): 7/8 up to 32 limbs, 3/4 from 48 limbs,
   where the high product saves more. */
#define MULHIGH_PAD_NUM 7
#define MULHIGH_PAD_DEN 8
#define MULHIGH_PAD_WIDE_N 48
#define MULHIGH_PAD_WIDE_NUM 3
#define MULHIGH_PAD_WIDE_DEN 4

FLINT_FORCE_INLINE int
_mulhigh_pad(mp_size_t m, mp_size_t n)
{
    return MULHIGH_PAD_DEN * m >= MULHIGH_PAD_NUM * n ||
        (n >= MULHIGH_PAD_WIDE_N && MULHIGH_PAD_WIDE_DEN * m >= MULHIGH_PAD_WIDE_NUM * n);
}

/* the full product computed on the stack (small) or in scratch, the top n
   limbs copied out; zn = xn + yn > n */
FLINT_FORCE_INLINE mp_limb_t
_mulhigh_via_full(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn, mp_size_t n)
{
    mp_size_t zn = xn + yn, d = zn - n - 1;
    mp_limb_t g;

    if (zn <= MULHIGH_STACK_LIMBS)
    {
        mp_limb_t t[MULHIGH_STACK_LIMBS];
        _mul_any(t, x, xn, y, yn);
        _flint_mpn_copy_small(res, t + d + 1, n);
        g = t[d];
    }
    else
    {
        mp_ptr t = flint_malloc(sizeof(mp_limb_t) * zn);
        _mul_any(t, x, xn, y, yn);
        flint_mpn_copyi(res, t + d + 1, n);
        g = t[d];
        flint_free(t);
    }
    return g;
}

/* scratch of the padded and peeled paths held on the stack up to this
   many output limbs (TMP_ALLOC and its size test cost about as much as a
   small kernel) */
#define MULHIGH_TRUNC_STACK_N 64

/* xn, yn <= n; res must not overlap x or y (inlined into the entry
   points, which have nothing else to do on the common path) */
FLINT_FORCE_INLINE mp_limb_t
_flint_mpn_mulhigh_trunc(mp_ptr res, mp_srcptr x, mp_size_t xn,
                         mp_srcptr y, mp_size_t yn, mp_size_t n)
{
    mp_size_t zn = xn + yn, d;
    mp_limb_t g, cy, stack[2 * MULHIGH_TRUNC_STACK_N];
    mp_ptr tmp;

    FLINT_ASSERT(xn >= 1 && yn >= 1 && xn <= n && yn <= n);

    if (xn == n && yn == n)
        return flint_mpn_mulhigh_n(res, x, y, n);

    /* the whole product fits in the output: exact */
    if (zn <= n)
    {
        _mul_any(res + n - zn, x, xn, y, yn);
        _flint_mpn_zero_small(res, n - zn);
        return 0;
    }

    d = zn - n - 1;

    /* small shapes, and shapes where nearly the whole product is wanted:
       compute all of it */
    if (n <= MULHIGH_FULL_PRODUCT_MAX_N || (zn <= MULHIGH_STACK_LIMBS && d <= 2))
        return _mulhigh_via_full(res, x, xn, y, yn, n);

#if FLINT_HAVE_FFT_SMALL
    if (FLINT_MAX(xn, yn) >= FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF && d >= 1)
    {
        if (FLINT_MIN(xn, yn) >= FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF)
        {
            /* the window of the product, its two lowest limbs diverted */
            mp_limb_t zlow[2];
            _mpn_ctx_mpn_mul_range_split(get_default_mpn_ctx(), res, zlow, 2,
                                         d - 1, zn, x, xn, y, yn);
            return zlow[1];
        }

        /* one operand below the FFT range: the full product, which the
           FFT computes faster than any split (measured) */
        return _mulhigh_via_full(res, x, xn, y, yn, n);
    }
#endif

    tmp = (n <= MULHIGH_TRUNC_STACK_N) ? stack : flint_malloc(sizeof(mp_limb_t) * 2 * n);

    if (_mulhigh_pad(FLINT_MIN(xn, yn), n))
    {
        /* nearly balanced: zero pad below to the balanced high product
           (same limbs of the product, same deficit) */
        mp_ptr pad = tmp;
        if (xn < n)
        {
            _flint_mpn_zero_small(pad, n - xn);
            _flint_mpn_copy_small(pad + n - xn, x, xn);
            x = pad;
            pad += n;
        }
        if (yn < n)
        {
            _flint_mpn_zero_small(pad, n - yn);
            _flint_mpn_copy_small(pad + n - yn, y, yn);
            y = pad;
        }
        g = flint_mpn_mulhigh_n(res, x, y, n);
    }
    else if (xn == n || yn == n)
    {
        /* one operand of full length: x = x_hi B^yn + x_lo, with x_hi * y
           landing exactly on the output */
        if (yn == n)
        {
            FLINT_SWAP(mp_srcptr, x, y);
            FLINT_SWAP(mp_size_t, xn, yn);
        }

        FLINT_ASSERT(xn == n && yn < n);

        _mul_any(res, x + yn, n - yn, y, yn);
        g = flint_mpn_mulhigh_n(tmp, x, y, yn);
        cy = mpn_add_n(res, res, tmp, yn);
        MPN_INCR_U(res + yn, n - yn, cy);
    }
    else
    {
        /* both shorter than the output and not nearly balanced: the full
           product (splitting off a high product does not pay in the
           Karatsuba range) */
        if (tmp != stack)
            flint_free(tmp);
        return _mulhigh_via_full(res, x, xn, y, yn, n);
    }

    if (tmp != stack)
        flint_free(tmp);
    return g;
}

/* x truncated to its top n limbs, the correction for the dropped limbs
   accumulated in (c1, c0) */
#define TRUNCATE_OPERAND(x, xn, ytop, n, c1, c0) \
    do { \
        if ((xn) > (n)) \
        { \
            mp_limb_t __h = n_mulhi((x)[(xn) - (n) - 1], (ytop)); \
            add_ssaaaa(c1, c0, c1, c0, UWORD(0), __h); \
            (x) += (xn) - (n); \
            (xn) = (n); \
        } \
    } while (0)

/* add the correction to the guard limb, carrying into the result */
#define APPLY_CORRECTION(res, n, g, c1, c0) \
    do { \
        if ((c1) | (c0)) \
        { \
            add_ssaaaa(c1, g, c1, g, UWORD(0), c0); \
            MPN_INCR_U(res, n, c1); \
        } \
    } while (0)

/*
    Aliasing: the balanced functions allow the output to be an operand;
    the general ones write the output while reading the operands on some
    paths, and after truncation the output buffer of an aliased operand
    overlaps its window without being equal to it, so an operand whose
    buffer is the output is copied (its window of at most n limbs: a
    negligible cost next to the product) before anything is written.
*/
mp_limb_t
flint_mpn_mulhigh(mp_ptr res, mp_srcptr x, mp_size_t xn,
                  mp_srcptr y, mp_size_t yn, mp_size_t n)
{
    mp_limb_t g, c1 = 0, c0 = 0;
    mp_srcptr xorig = x, yorig = y;
    mp_ptr tmp = NULL;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(xn >= 1 && yn >= 1);

    if (xn == n && yn == n)
        return flint_mpn_mulhigh_n(res, x, y, n);

    /* no truncation, no aliasing: nothing to set up */
    if (xn <= n && yn <= n && res != x && res != y)
        return _flint_mpn_mulhigh_trunc(res, x, xn, y, yn, n);

    TRUNCATE_OPERAND(x, xn, y[yn - 1], n, c1, c0);
    TRUNCATE_OPERAND(y, yn, x[xn - 1], n, c1, c0);

    TMP_START;
    if (res == xorig || res == yorig)
    {
        /* the windows of the aliased operands, copied */
        mp_ptr t = tmp = TMP_ALLOC(sizeof(mp_limb_t) * (xn + yn));
        if (res == xorig)
        {
            flint_mpn_copyi(t, x, xn);
            x = t;
            t += xn;
        }
        if (res == yorig)
        {
            if (yorig == xorig && yn == xn)
                y = x;
            else
            {
                flint_mpn_copyi(t, y, yn);
                y = t;
            }
        }
    }

    if (xn == n && yn == n)
        g = flint_mpn_mulhigh_n(res, x, y, n);
    else
        g = _flint_mpn_mulhigh_trunc(res, x, xn, y, yn, n);

    TMP_END;
    APPLY_CORRECTION(res, n, g, c1, c0);
    return g;
}

/* xn <= n; res must not overlap x */
FLINT_FORCE_INLINE mp_limb_t
_flint_mpn_sqrhigh_trunc(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_size_t zn = 2 * xn, d;
    mp_limb_t g, stack[MULHIGH_TRUNC_STACK_N];
    mp_ptr tmp;

    FLINT_ASSERT(xn >= 1 && xn <= n);

    if (xn == n)
        return flint_mpn_sqrhigh_n(res, x, n);

    if (zn <= n)
    {
        flint_mpn_sqr(res + n - zn, x, xn);
        _flint_mpn_zero_small(res, n - zn);
        return 0;
    }

    d = zn - n - 1;

    if (n <= MULHIGH_FULL_PRODUCT_MAX_N || (zn <= MULHIGH_STACK_LIMBS && d <= 2))
        return _mulhigh_via_full(res, x, xn, x, xn, n);

#if FLINT_HAVE_FFT_SMALL
    if (xn >= FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF && d >= 1)
    {
        mp_limb_t zlow[2];
        _mpn_ctx_mpn_mul_range_split(get_default_mpn_ctx(), res, zlow, 2,
                                     d - 1, zn, x, xn, x, xn);
        return zlow[1];
    }
#endif

    if (!_mulhigh_pad(xn, n))
        return _mulhigh_via_full(res, x, xn, x, xn, n);

    tmp = (n <= MULHIGH_TRUNC_STACK_N) ? stack : flint_malloc(sizeof(mp_limb_t) * n);
    _flint_mpn_zero_small(tmp, n - xn);
    _flint_mpn_copy_small(tmp + n - xn, x, xn);
    g = flint_mpn_sqrhigh_n(res, tmp, n);
    if (tmp != stack)
        flint_free(tmp);
    return g;
}

mp_limb_t
flint_mpn_sqrhigh(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_limb_t g, c1 = 0, c0 = 0;
    mp_srcptr xorig = x;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(xn >= 1);

    if (xn == n)
        return flint_mpn_sqrhigh_n(res, x, n);

    if (xn < n && res != x)
        return _flint_mpn_sqrhigh_trunc(res, x, xn, n);

    if (xn > n)
    {
        /* both cross terms: 2 * mulhi(x_{-1}, x_top) */
        mp_limb_t h = n_mulhi(x[xn - n - 1], x[xn - 1]);
        add_ssaaaa(c1, c0, c1, c0, UWORD(0), h);
        add_ssaaaa(c1, c0, c1, c0, UWORD(0), h);
        x += xn - n;
        xn = n;
    }

    TMP_START;
    if (res == xorig)
    {
        mp_ptr t = TMP_ALLOC(sizeof(mp_limb_t) * xn);
        flint_mpn_copyi(t, x, xn);
        x = t;
    }

    if (xn == n)
        g = flint_mpn_sqrhigh_n(res, x, n);
    else
        g = _flint_mpn_sqrhigh_trunc(res, x, xn, n);

    TMP_END;
    APPLY_CORRECTION(res, n, g, c1, c0);
    return g;
}

/* The padded normalised kernel for an operand short by one limb (two for
   larger n) beats the full product with a shift from this many limbs
   (measured on Zen 3: equal at 8, kernel ahead at 9). */
#define MULHIGH_NORMALISED_PAD_MIN_N 8

/* add the truncation correction (c1, c0) at the guard limb d = zn - n - 1
   of the full product t of the truncated operands (no overflow: the sum
   is at most the product of the untruncated operands) */
static void
_add_correction_to_product(mp_ptr t, mp_size_t zn, mp_size_t n, mp_limb_t c1, mp_limb_t c0)
{
    mp_size_t d = zn - n - 1;

    if (c1 | c0)
    {
        mpn_add_1(t + d, t + d, zn - d, c0);
        if (c1)
            mpn_add_1(t + d + 1, t + d + 1, zn - d - 1, c1);
    }
}

/* The top n limbs of the product t (zn limbs, zn > n), normalised: the
   copy and the shift in one pass. */
static mp_limb_pair_t
_copy_high_normalised(mp_ptr res, mp_srcptr t, mp_size_t zn, mp_size_t n)
{
    mp_size_t d = zn - n - 1, i;
    mp_limb_pair_t ret;

    if (t[zn - 1] >> (FLINT_BITS - 1))
    {
        flint_mpn_copyi(res, t + d + 1, n);
        ret.m1 = t[d];
        ret.m2 = 0;
    }
    else
    {
        for (i = 0; i < n; i++)
            res[i] = (t[d + 1 + i] << 1) | (t[d + i] >> (FLINT_BITS - 1));
        ret.m1 = (t[d] << 1) | ((d > 0) ? (t[d - 1] >> (FLINT_BITS - 1)) : 0);
        ret.m2 = 1;
    }

    return ret;
}

static mp_limb_pair_t
_normalise(mp_ptr res, mp_size_t n, mp_limb_t g)
{
    mp_limb_pair_t ret;

    ret.m1 = g;

    if (res[n - 1] >> (FLINT_BITS - 1))
    {
        ret.m2 = 0;
    }
    else
    {
        ret.m2 = 1;
        mpn_lshift(res, res, n, 1);
        res[0] |= (g >> (FLINT_BITS - 1));
        ret.m1 = g << 1;
    }

    return ret;
}

/* Add the truncation correction (c1, c0), in units of the unshifted guard
   limb, to a result already normalised in registers. If the kernel shifted,
   the exact value may nevertheless have had its top bit set (the deficit
   hid it), in which case the corrected shifted value overflows the n limbs:
   then shift back. */
static mp_limb_pair_t
_apply_correction_normalised(mp_ptr res, mp_size_t n, mp_limb_pair_t ret,
                             mp_limb_t c1, mp_limb_t c0)
{
    mp_limb_t g = ret.m1, cy;

    if (ret.m2)
    {
        c1 = (c1 << 1) | (c0 >> (FLINT_BITS - 1));
        c0 <<= 1;
    }

    add_ssaaaa(c1, g, c1, g, UWORD(0), c0);
    cy = mpn_add_1(res, res, n, c1);

    if (cy)
    {
        FLINT_ASSERT(ret.m2);
        g = (res[0] << (FLINT_BITS - 1)) | (g >> 1);
        mpn_rshift(res, res, n, 1);
        res[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
        ret.m2 = 0;
    }

    ret.m1 = g;
    return ret;
}

mp_limb_pair_t
flint_mpn_mulhigh_normalised(mp_ptr res, mp_srcptr x, mp_size_t xn,
                             mp_srcptr y, mp_size_t yn, mp_size_t n)
{
    mp_limb_t g, c1 = 0, c0 = 0;
    mp_limb_pair_t ret;

    mp_srcptr xorig = x, yorig = y;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(xn >= 1 && yn >= 1);

    if (xn == n && yn == n)
        return flint_mpn_mulhigh_n_normalised(res, x, y, n);

    TRUNCATE_OPERAND(x, xn, y[yn - 1], n, c1, c0);
    TRUNCATE_OPERAND(y, yn, x[xn - 1], n, c1, c0);

    TMP_START;
    if (res == xorig || res == yorig)
    {
        mp_ptr t = TMP_ALLOC(sizeof(mp_limb_t) * (xn + yn));
        if (res == xorig)
        {
            flint_mpn_copyi(t, x, xn);
            x = t;
            t += xn;
        }
        if (res == yorig)
        {
            if (yorig == xorig && yn == xn)
                y = x;
            else
            {
                flint_mpn_copyi(t, y, yn);
                y = t;
            }
        }
    }

    if (xn == n && yn == n)
    {
        ret = flint_mpn_mulhigh_n_normalised(res, x, y, n);
    }
#if FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH > 0
    else if (n <= FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH && n >= MULHIGH_NORMALISED_PAD_MIN_N
        && FLINT_MAX(xn, yn) == n && FLINT_MIN(xn, yn) + 1 + (n > MULHIGH_NORMALISED_PAD_MIN_N) >= n)
    {
        /* one operand of full length, the other short by one or two limbs:
           zero padding it for the normalised kernel (which shifts in
           registers) beats the full product followed by a shift */
        mp_limb_t b[FLINT_MPN_MULHIGH_NORMALISED_FUNC_TAB_WIDTH];
        if (xn < n)
        {
            _flint_mpn_zero_small(b, n - xn);
            _flint_mpn_copy_small(b + n - xn, x, xn);
            ret = flint_mpn_mulhigh_normalised_func_tab[n](res, b, y);
        }
        else
        {
            _flint_mpn_zero_small(b, n - yn);
            _flint_mpn_copy_small(b + n - yn, y, yn);
            ret = flint_mpn_mulhigh_normalised_func_tab[n](res, x, b);
        }
    }
#endif
    else if (n <= MULHIGH_FULL_PRODUCT_MAX_N && xn + yn > n)
    {
        /* small: the full product, copied with the shift */
        mp_limb_t t[MULHIGH_STACK_LIMBS];
        _mul_any(t, x, xn, y, yn);
        _add_correction_to_product(t, xn + yn, n, c1, c0);
        TMP_END;
        return _copy_high_normalised(res, t, xn + yn, n);
    }
    else
    {
        g = _flint_mpn_mulhigh_trunc(res, x, xn, y, yn, n);
        TMP_END;
        APPLY_CORRECTION(res, n, g, c1, c0);
        return _normalise(res, n, g);
    }

    TMP_END;

    if (c1 | c0)
        ret = _apply_correction_normalised(res, n, ret, c1, c0);

    return ret;
}

mp_limb_pair_t
flint_mpn_sqrhigh_normalised(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
{
    mp_limb_t g, c1 = 0, c0 = 0;
    mp_limb_pair_t ret;

    mp_srcptr xorig = x;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(xn >= 1);

    if (xn == n)
        return flint_mpn_sqrhigh_n_normalised(res, x, n);

    if (xn > n)
    {
        mp_limb_t h = n_mulhi(x[xn - n - 1], x[xn - 1]);
        add_ssaaaa(c1, c0, c1, c0, UWORD(0), h);
        add_ssaaaa(c1, c0, c1, c0, UWORD(0), h);
        x += xn - n;
        xn = n;
    }

    TMP_START;
    if (res == xorig)
    {
        mp_ptr t = TMP_ALLOC(sizeof(mp_limb_t) * xn);
        flint_mpn_copyi(t, x, xn);
        x = t;
    }

    if (xn == n)
    {
        ret = flint_mpn_sqrhigh_n_normalised(res, x, n);
    }
#if FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH > 0
    else if (n <= FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH && n >= MULHIGH_NORMALISED_PAD_MIN_N
        && xn + 1 + (n > MULHIGH_NORMALISED_PAD_MIN_N) >= n)
    {
        mp_limb_t b[FLINT_MPN_SQRHIGH_NORMALISED_FUNC_TAB_WIDTH];
        _flint_mpn_zero_small(b, n - xn);
        _flint_mpn_copy_small(b + n - xn, x, xn);
        ret = flint_mpn_sqrhigh_normalised_func_tab[n](res, b);
    }
#endif
    else if (n <= MULHIGH_FULL_PRODUCT_MAX_N && 2 * xn > n)
    {
        mp_limb_t t[MULHIGH_STACK_LIMBS];
        flint_mpn_sqr(t, x, xn);
        _add_correction_to_product(t, 2 * xn, n, c1, c0);
        TMP_END;
        return _copy_high_normalised(res, t, 2 * xn, n);
    }
    else
    {
        g = _flint_mpn_sqrhigh_trunc(res, x, xn, n);
        TMP_END;
        APPLY_CORRECTION(res, n, g, c1, c0);
        return _normalise(res, n, g);
    }

    TMP_END;

    if (c1 | c0)
        ret = _apply_correction_normalised(res, n, ret, c1, c0);

    return ret;
}

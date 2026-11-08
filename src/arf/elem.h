/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal definitions for the elementary functions of arf, shared with
   nfloat (which includes this header). */

#ifndef ARF_ELEM_H
#define ARF_ELEM_H

#include <math.h>
#include "longlong.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "arf.h"
#include "../mpn_extras/fp_small.h"
#include "../mp_real/impl.h"

/*
    Conventions (mp_real.h): the evaluation is done by the _fp functions
    of mp_real on fixed-point frames at nk kernel limbs, with error
    counts in ulps of the result.  The result is then rounded: nfloat
    truncates it in its default rounding mode and rounds outward in its
    directed modes (_nfloat_set_mpn_err); arf rounds it in the relaxed
    modes ARF_RND_FAST and ARF_RND_ACCURATE, and in the strict modes
    certifies the rounding from the bits between the rounding point and
    the error bound, recomputing at a higher precision when they do not
    decide it (_arf_set_mpn_err, and the loops of the arf functions).

    The target is a relative accuracy of 2^-p, p the precision of the
    arf call; the kernels run at enough limbs for p plus
    ARF_ELEM_GUARD_BITS bits in the relaxed modes (the kernels' errors,
    below 2^13 ulps, must stay below an eighth of an ulp for the 1 ulp
    bound of ARF_RND_FAST and below 2^-7 ulp for the 0.51 ulp bound of
    ARF_RND_ACCURATE: _arf_set_mpn_err checks and asks for a retry
    otherwise), and ARF_ELEM_STRICT_GUARD_BITS in the strict modes, which
    leaves a window of bits for the certification.
*/

#define ARF_ELEM_GUARD_BITS 20
#define ARF_ELEM_STRICT_GUARD_BITS 24

/* the return value of a fast path or a certification that does not
   decide the result */
#define ARF_ELEM_UNCERTAIN (-1)

/* the arf functions compute exp(x) resp. the trigonometric functions of
   x for |x| up to 2^ARF_EXP_MAGLIM(prec) resp. 2^ARF_TRIG_MAGLIM(prec)
   (the limits of arb) and return nan beyond */
#define ARF_EXP_MAGLIM(prec) FLINT_MAX(128, 2 * (prec))
#define ARF_TRIG_MAGLIM(prec) MP_REAL_TRIG_MAGLIM(prec)

FLINT_FORCE_INLINE int
_arf_rnd_is_strict(arf_rnd_t rnd)
{
    return rnd < ARF_RND_FAST;
}

/* the bits for the kernels (mp_real.h) for an absolute accuracy of
   2^-(prec + extra) in the mode rnd (the guard bits above), and the
   fixed-point limbs for them */
FLINT_FORCE_INLINE slong
_arf_elem_bits(slong prec, arf_rnd_t rnd, slong extra)
{
    slong guard = _arf_rnd_is_strict(rnd) ? ARF_ELEM_STRICT_GUARD_BITS : ARF_ELEM_GUARD_BITS;
    return prec + guard + extra;
}

FLINT_FORCE_INLINE slong
_arf_elem_limbs(slong prec, arf_rnd_t rnd, slong extra)
{
    return _mp_real_elem_limbs(_arf_elem_bits(prec, rnd, extra));
}

/* Result assembly for arf **************************************************/

/* whether the bits [lo, hi) of y are all zero (ones = 0) resp. all one
   (ones = 1); hi > lo, and y has at least hi bits */
FLINT_FORCE_INLINE int
_arf_mpn_bits_uniform(nn_srcptr y, slong lo, slong hi, int ones)
{
    slong i, qlo = lo / FLINT_BITS, qhi = (hi - 1) / FLINT_BITS;
    ulong v, mask;

    for (i = qlo; i <= qhi; i++)
    {
        mask = ~UWORD(0);
        if (i == qlo && lo % FLINT_BITS != 0)
            mask &= ~UWORD(0) << (lo % FLINT_BITS);
        if (i == qhi && hi % FLINT_BITS != 0)
            mask &= ~UWORD(0) >> (FLINT_BITS - hi % FLINT_BITS);
        v = y[i] & mask;
        if (v != (ones ? mask : 0))
            return 0;
    }

    return 1;
}

/* res = (-1)^sgnbit V for V = (hi, lo) 2^(exp - 2 FLINT_BITS) + (bits
   below: half the first, sticky whether any other is set), hi with its
   top bit set, rounded to prec <= 2 FLINT_BITS bits in the sign-free
   mode rnd (ARF_RND_DOWN, ARF_RND_UP or ARF_RND_NEAR); the exponent exp
   must be a small fmpz. */
FLINT_FORCE_INLINE void
_arf_set_frame2(arf_t res, ulong hi, ulong lo, int half, int sticky, slong exp, int sgnbit, slong prec, arf_rnd_t rnd)
{
    slong sb = 2 * FLINT_BITS - prec;
    ulong h, st, lsb;
    int up;

    if (sb == 0)
    {
        h = half;
        st = sticky;
        lsb = lo & 1;
    }
    else if (sb <= FLINT_BITS)
    {
        h = (lo >> (sb - 1)) & 1;
        st = sticky | half | ((sb > 1) && (lo << (FLINT_BITS + 1 - sb)) != 0);
        lo = (sb == FLINT_BITS) ? 0 : (lo & (~UWORD(0) << sb));
        lsb = (sb == FLINT_BITS) ? (hi & 1) : ((lo >> sb) & 1);
    }
    else
    {
        slong s2 = sb - FLINT_BITS;
        h = (hi >> (s2 - 1)) & 1;
        st = sticky | half | (lo != 0) | ((s2 > 1) && (hi << (FLINT_BITS + 1 - s2)) != 0);
        hi &= ~UWORD(0) << s2;
        lo = 0;
        lsb = (hi >> s2) & 1;
    }

    if (rnd == ARF_RND_DOWN)
        up = 0;
    else if (rnd == ARF_RND_UP)
        up = (h | st) != 0;
    else
        up = h && (st || lsb);

    if (up)
    {
        if (sb < FLINT_BITS)
            add_ssaaaa(hi, lo, hi, lo, UWORD(0), UWORD(1) << sb);
        else
            hi += UWORD(1) << (sb - FLINT_BITS);

        if (hi == 0)
        {
            /* a carry out: the next power of two */
            hi = LIMB_TOP;
            lo = 0;
            exp++;
        }
    }

    ARF_DEMOTE(res);
    if (lo == 0)
    {
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(res)[0] = hi;
    }
    else
    {
        ARF_XSIZE(res) = ARF_MAKE_XSIZE(2, sgnbit);
        ARF_NOPTR_D(res)[0] = lo;
        ARF_NOPTR_D(res)[1] = hi;
    }
    _fmpz_demote(ARF_EXPREF(res));
    ARF_EXP(res) = exp;
}

/* Sets res to (-1)^sgnbit V where V is a real number with
   (Y - 2^eblo) 2^e < V < (Y + 2^ebhi) 2^e for the integer Y = (y, len) > 0,
   rounded to prec bits in the mode rnd (y needs room for len + 1 limbs
   and is destroyed; the exponent e is a plain word, so the callers handle
   any exponent of V beyond the word range).  A negative eblo resp. ebhi
   (with the other one nonnegative) means that V > Y 2^e resp.
   V < Y 2^e; both negative mean V = Y 2^e exactly, which is then rounded
   exactly (and its inexact flag returned).  The relaxed modes round
   Y + 2^ebhi truncated (ARF_RND_FAST: within an ulp of the value when
   the error bounds stay below an eighth of an ulp) resp. Y to nearest
   (ARF_RND_ACCURATE: within 0.51 ulp when they stay below 2^-7 ulp),
   and return ARF_ELEM_UNCERTAIN when the bounds are larger.
   The strict modes need the bits of Y between the error and the
   rounding point to decide the rounding of V (elem.c explains the
   rules); V is assumed inexact at prec bits when an error is given.
   Returns 1 (the inexact flag), or ARF_ELEM_UNCERTAIN without writing
   res when the window does not decide the rounding (in any mode when
   Y = 0, as an approximation of a nonzero value). */
int _arf_set_mpn_err(arf_t res, nn_ptr y, slong len, slong e, slong eblo, slong ebhi, int sgnbit, slong prec, arf_rnd_t rnd);

/* The inverse trigonometric and hyperbolic functions at one and two
   limbs by in-register compositions (arf/elem_arc_fast.c), for the
   mantissa (d, n) (n = 1 or 2; a shorter mantissa padded with a zero
   low limb) with exponent e, normal x in the domain (|x| != 1 where that
   is special): the result (-1)^sgnbit (T, w + 1) 2^(k - FLINT_BITS w),
   w = n + 1, with the integral limb T[w], within err units; returns
   ARF_ELEM_UNCERTAIN when the path does not apply (an intermediate
   losing relative accuracy, |x| >= 2^62, |x| < 2^-31 resp. 2^-39). */
#define ARF_ARC_ASIN 0
#define ARF_ARC_ACOS 1
#define ARF_ARC_ATAN 2
#define ARF_ARC_ACOT 3
#define ARF_ARC_ASEC 4
#define ARF_ARC_ACSC 5
#define ARF_ARC_ASINH 6
#define ARF_ARC_ACOSH 7
#define ARF_ARC_ATANH 8
#define ARF_ARC_ACOTH 9
#define ARF_ARC_ASECH 10
#define ARF_ARC_ACSCH 11
typedef struct
{
    ulong T[4];
    slong k;
    ulong err;
    int sgnbit;
}
_arf_arc_out_struct;
#if FLINT_BITS == 64
int _arf_arc_fast(_arf_arc_out_struct * res, nn_srcptr d, slong n, slong e, int xsgn, int f, int div_pi);
#endif


/* res = (-1)^sgnbit Q for the quotient Q = V / D (invert = 0) or D / V
   (invert = 1) of V = (y, len) 2^ue, known within err units, by the
   exact D = (d, dn) 2^(de - FLINT_BITS dn) (a normalised mantissa; D = 1
   for d = NULL): V normalised to wn limbs (the working precision of the
   caller's retry, so that the retries make progress) and one approximate
   division (_arf_mpn_div_approx), then rounded as _arf_set_mpn_err does;
   ARF_ELEM_UNCERTAIN for V = 0 or an error count the normalisation
   cannot scale */
int _arf_set_mpn_err_quot(arf_t res, nn_srcptr y, slong len, slong ue, ulong err, nn_srcptr d, slong dn, slong de, int invert, int sgnbit, slong wn, slong prec, arf_rnd_t rnd);

/* res = (-1)^rsgn exp(x) 2^shift (within elo units below and ehi above
   of the kernel's frame, on top of its error) for x = (-1)^sgnbit (d, xn)
   2^(e - FLINT_BITS xn) with MP_REAL_EXP_FP_MAX_EXP < e (arf/exp.c), by
   _mp_real_exp_fp_huge with the bignum quotient as the exponent of the
   result; nan beyond e = ARF_EXP_MAGLIM(prec) */
int _arf_exp_huge(arf_t res, nn_srcptr d, slong xn, slong e, int sgnbit, ulong elo, ulong ehi, slong shift, int rsgn, slong prec, arf_rnd_t rnd);

/* res = (1/x)(1 + delta) with -2^-P < delta < 0 (side = -1) resp.
   0 < delta < 2^-P (side = 1), P >= prec + 3: the reciprocal of the
   mantissa floored exactly (the remainder tells whether 1/x is a dyadic
   number), or ARF_ELEM_UNCERTAIN when the bits of 1/x between the
   rounding point and 2^-P do not decide the rounding (which they do for
   P >= prec + FLINT_BITS (xn + 2), the expansion of 1/x having no run of
   equal bits that long) */
int _arf_set_tiny_recip_rel_err(arf_t res, const arf_t x, int side, slong P, slong prec, arf_rnd_t rnd);

/* the bound 2^eb > err of an error count in ulps (eb = -1 for err = 0,
   an exact value) */
#define ARF_ELEM_EB(err) (((err) == 0) ? WORD(-1) : (slong) FLINT_BIT_COUNT(err))

/* res = (-1)^sgnbit (1 + delta) with 0 < delta < 2^(-prec - 1) (dsgn = 0)
   resp. -2^(-prec - 1) < delta < 0 (dsgn = 1), correctly rounded to prec
   bits in the mode rnd (1 or its neighbour); returns 1 */
int _arf_one_plus_tiny(arf_t res, int sgnbit, int dsgn, slong prec, arf_rnd_t rnd);

/* res = x (1 + delta) with -2^-P < delta < 0 (side = -1) resp.
   0 < delta < 2^-P (side = 1), P >= prec + 3, correctly rounded to prec
   bits in the mode rnd for the (nonzero, finite) arf x; returns 1, or
   ARF_ELEM_UNCERTAIN when x has more than P bits and its bits between
   the rounding point and 2^-P do not decide the rounding */
int _arf_set_tiny_rel_err(arf_t res, const arf_t x, int side, slong P, slong prec, arf_rnd_t rnd);

/* the short series of a small exact argument |x| = (d, xn)
   2^(e - FLINT_BITS xn) (arf/elem.c): x + c x^3 for the odd functions,
   1 + c x^2 for the even ones, x -+ x^2/2 + x^3/3 for log1p (the tail
   below 2^(4e - 1) for |x| < 1/4), in fixed point at n limbs (nk + 2 for
   the kernels' nk), Y (n + 2 limbs) receiving n + 1 limbs with the unit
   exponent ue and the error counts elo, ehi; the caller adds the tail
   (the shift -2e is at least 2); returns whether the function is even
   (no sign carried) */
#define ARF_SMALL_SERIES_SIN 0
#define ARF_SMALL_SERIES_TAN 1
#define ARF_SMALL_SERIES_ATAN 2
#define ARF_SMALL_SERIES_SINH 3
#define ARF_SMALL_SERIES_ASINH 4
#define ARF_SMALL_SERIES_ATANH 5
#define ARF_SMALL_SERIES_COS 6
#define ARF_SMALL_SERIES_COSH 7
#define ARF_SMALL_SERIES_LOG1P_POS 8    /* |log(1 + x)| = x - x^2/2 + x^3/3 for x > 0 */
#define ARF_SMALL_SERIES_LOG1P_NEG 9    /* |log(1 + x)| = |x| + x^2/2 + |x|^3/3 for x < 0 */
int _arf_small_arg_series_frame(nn_ptr Y, slong * ue, ulong * elo, ulong * ehi, nn_srcptr d, slong xn, slong e, int kind, slong n);

/* the series rounded into res for x = (-1)^sgnbit (d, xn) 2^(e - FLINT_BITS xn),
   the tail below 2^(5e + tail) (odd) resp. 2^(4e + tail) (even), with
   the rules of _arf_set_mpn_err (ARF_ELEM_UNCERTAIN when the window
   does not decide) */
int _arf_small_arg_series(arf_t res, nn_srcptr d, slong xn, slong e, int sgnbit, int kind, slong tail, slong prec, arf_rnd_t rnd);

/* the leading zero bits of |x - 1| for x = M 2^e, e in {0, 1}, or WORD_MAX
   for x = 1 */
FLINT_FORCE_INLINE slong
_arf_log_lzb(nn_srcptr d, slong xn, slong e)
{
    slong j, z;

    if (e == 1)
    {
        /* the fraction of 2M */
        ulong top = d[xn - 1] << 1;

        if (top != 0)
            return flint_clz(top);
        for (j = xn - 2; j >= 0 && d[j] == 0; j--)
            ;
        if (j < 0)
            return WORD_MAX;
        /* d[j] contributes its top bit to the limb above */
        z = FLINT_BITS * (xn - 1 - j) + flint_clz(d[j]) - 1;
        return z;
    }
    else
    {
        /* 1 - M: the leading ones of M */
        z = 0;
        for (j = xn - 1; j >= 0 && d[j] == UWORD_MAX; j--)
            z += FLINT_BITS;
        if (j >= 0)
            z += flint_clz(~d[j]);
        return z - 1;
    }
}

/* the functions (which) of the hyperbolic fast paths, also named by the
   callers in arb, nfloat and acf on 32-bit machines, where the paths
   below (64-bit only) are not compiled */
#define ARF_SMALL_HYP_TANH 0
#define ARF_SMALL_HYP_COTH 1
#define ARF_SMALL_HYP_SECH 2
#define ARF_SMALL_HYP_CSCH 3

#if FLINT_BITS == 64
/* The finishing step of the in-register fast paths at one and two limbs
   (64-bit only): res = (-1)^sgnbit Y 2^(k - FLINT_BITS (len - 1)) for
   Y = (y2, y1, y0) resp. (y3, y2, y1, y0) with the integral limb y2 resp.
   y3 below 2^(FLINT_BITS - 1), the top two limbs not both zero, known
   within err units: rounded to prec <= FLINT_BITS resp. 2 FLINT_BITS bits
   in the mode rnd by the rules of _arf_set_mpn_err (two-sided error), or
   ARF_ELEM_UNCERTAIN.  The exponent k + ... must stay within the word
   range. */
/* the bits [lo, lo + w) of the frame (f2, f1, f0), w <= FLINT_BITS */
FLINT_FORCE_INLINE ulong
_frame3_bits(ulong f2, ulong f1, ulong f0, slong lo, slong w)
{
    ulong v;

    if (lo >= 2 * FLINT_BITS)
        v = f2 >> (lo - 2 * FLINT_BITS);
    else if (lo == FLINT_BITS)
        v = f1;
    else if (lo > FLINT_BITS)
        v = (f1 >> (lo - FLINT_BITS)) | (f2 << (2 * FLINT_BITS - lo));
    else if (lo == 0)
        v = f0;
    else
        v = (f0 >> lo) | (f1 << (FLINT_BITS - lo));

    return (w == FLINT_BITS) ? v : (v & ((UWORD(1) << w) - 1));
}

FLINT_FORCE_INLINE int
_arf_small_set_1(arf_t res, ulong y2, ulong y1, ulong y0, slong k, ulong err, int sgnbit, slong prec, arf_rnd_t rnd)
{
    ulong hi, lo, rest;
    unsigned int c;
    slong exp;

    /* ARF_RND_FAST truncates the upper bound */
    if (rnd == ARF_RND_FAST)
        add_sssaaaaaa(y2, y1, y0, y2, y1, y0, UWORD(0), UWORD(0), err);

    /* normalise to the 192-bit frame (hi, lo, rest); the error scales
       by 2^c */
    if (y2 != 0)
    {
        c = flint_clz(y2);
        if (c == 0)
        {
            hi = y2; lo = y1; rest = y0;
        }
        else
        {
            hi = (y2 << c) | (y1 >> (FLINT_BITS - c));
            lo = (y1 << c) | (y0 >> (FLINT_BITS - c));
            rest = y0 << c;
        }
        exp = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y1);
        if (c == 0)
        {
            hi = y1; lo = y0;
        }
        else
        {
            hi = (y1 << c) | (y0 >> (FLINT_BITS - c));
            lo = y0 << c;
        }
        rest = 0;
        exp = k - (slong) c;
        /* (the error scales by the whole shift of the frame) */
        c += FLINT_BITS;
    }

    if (_arf_rnd_is_strict(rnd))
    {
        /* the window [eb, top) below the half bit at top = 191 - prec,
           at most 64 bits of it from the top (the fast paths have at
           least 50), neither all zero nor all one */
        slong eb = (slong) FLINT_BIT_COUNT(err) + c;
        slong top = 3 * FLINT_BITS - 1 - prec;
        slong w = top - eb;
        ulong v, mask;

        if (w < 1)
            return ARF_ELEM_UNCERTAIN;
        w = FLINT_MIN(w, FLINT_BITS);
        v = _frame3_bits(hi, lo, rest, top - w, w);
        mask = (w == FLINT_BITS) ? ~UWORD(0) : ((UWORD(1) << w) - 1);
        if (v == 0 || v == mask)
            return ARF_ELEM_UNCERTAIN;

        /* decided: rounds like the approximation with a sticky bit */
        _arf_set_frame2(res, hi, lo, rest >> (FLINT_BITS - 1), 1, exp, sgnbit, prec,
            _arf_rnd_resolve(rnd, sgnbit));
    }
    else
    {
        /* the error (in the frame, with the added bound of FAST) below an
           eighth resp. 2^-7 of an ulp */
        slong eb = (slong) FLINT_BIT_COUNT(err) + c + 1;
        slong top = 3 * FLINT_BITS - 1 - prec;

        if (eb > top - ((rnd == ARF_RND_FAST) ? 3 : 7))
            return ARF_ELEM_UNCERTAIN;
        _arf_set_frame2(res, hi, lo, rest >> (FLINT_BITS - 1), (rest << 1) != 0, exp, sgnbit, prec,
            (rnd == ARF_RND_FAST) ? ARF_RND_DOWN : ARF_RND_NEAR);
    }

    return 1;
}

FLINT_FORCE_INLINE int
_arf_small_set_2(arf_t res, ulong y3, ulong y2, ulong y1, ulong y0, slong k, ulong err, int sgnbit, slong prec, arf_rnd_t rnd)
{
    ulong t[4];
    unsigned int c;
    slong exp;

    if (rnd == ARF_RND_FAST)
        add_ssssaaaaaaaa(y3, y2, y1, y0, y3, y2, y1, y0, UWORD(0), UWORD(0), UWORD(0), err);

    /* the 256-bit frame t, normalised */
    if (y3 != 0)
    {
        c = flint_clz(y3);
        if (c == 0)
        {
            t[3] = y3; t[2] = y2; t[1] = y1; t[0] = y0;
        }
        else
        {
            t[3] = (y3 << c) | (y2 >> (FLINT_BITS - c));
            t[2] = (y2 << c) | (y1 >> (FLINT_BITS - c));
            t[1] = (y1 << c) | (y0 >> (FLINT_BITS - c));
            t[0] = y0 << c;
        }
        exp = k + FLINT_BITS - c;
    }
    else
    {
        c = flint_clz(y2);
        if (c == 0)
        {
            t[3] = y2; t[2] = y1; t[1] = y0;
        }
        else
        {
            t[3] = (y2 << c) | (y1 >> (FLINT_BITS - c));
            t[2] = (y1 << c) | (y0 >> (FLINT_BITS - c));
            t[1] = y0 << c;
        }
        t[0] = 0;
        exp = k - (slong) c;
        c += FLINT_BITS;
    }

    if (_arf_rnd_is_strict(rnd))
    {
        /* the window [eb, top) below the half bit at top = 255 - prec */
        slong eb = (slong) FLINT_BIT_COUNT(err) + c;
        slong top = 4 * FLINT_BITS - 1 - prec;

        if (top - eb < 1 || _arf_mpn_bits_uniform(t, eb, top, 0) || _arf_mpn_bits_uniform(t, eb, top, 1))
            return ARF_ELEM_UNCERTAIN;

        _arf_set_frame2(res, t[3], t[2], t[1] >> (FLINT_BITS - 1), 1, exp, sgnbit, prec,
            _arf_rnd_resolve(rnd, sgnbit));
    }
    else
    {
        slong eb = (slong) FLINT_BIT_COUNT(err) + c + 1;
        slong top = 4 * FLINT_BITS - 1 - prec;

        if (eb > top - ((rnd == ARF_RND_FAST) ? 3 : 7))
            return ARF_ELEM_UNCERTAIN;
        _arf_set_frame2(res, t[3], t[2], t[1] >> (FLINT_BITS - 1), ((t[1] << 1) | t[0]) != 0, exp, sgnbit, prec,
            (rnd == ARF_RND_FAST) ? ARF_RND_DOWN : ARF_RND_NEAR);
    }

    return 1;
}

/* Fixed-point helpers of the fast paths (64-bit). */

/* |x| = d 2^(e - FLINT_BITS) as (xi, x1, x0) with the binary point after
   the integral limb xi, exactly, for -FLINT_BITS <= e < FLINT_BITS */
FLINT_FORCE_INLINE void
_arf_small_fixed_1(ulong * xi, ulong * x1, ulong * x0, ulong d, slong e)
{
    if (e > 0)
    {
        *xi = d >> (FLINT_BITS - e);
        *x1 = d << e;
        *x0 = 0;
    }
    else if (e == 0)
    {
        *xi = 0; *x1 = d; *x0 = 0;
    }
    else if (e > -FLINT_BITS)
    {
        *xi = 0;
        *x1 = d >> (-e);
        *x0 = d << (FLINT_BITS + e);
    }
    else
    {
        *xi = 0; *x1 = 0; *x0 = d;
    }
}

/* |x| = (d1, d0) 2^(e - 2 FLINT_BITS) as (xi, x2, x1, x0), binary point
   after xi, for -2 FLINT_BITS <= e < FLINT_BITS: exact for e >=
   -FLINT_BITS, truncated below B^-3 otherwise */
FLINT_FORCE_INLINE void
_arf_small_fixed_2(ulong * xi, ulong * x2, ulong * x1, ulong * x0, ulong d1, ulong d0, slong e)
{
    unsigned int s;

    *xi = *x2 = *x1 = *x0 = 0;
    if (e > 0)
    {
        *xi = d1 >> (FLINT_BITS - e);
        *x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        *x1 = d0 << e;
    }
    else if (e == 0)
    {
        *x2 = d1; *x1 = d0;
    }
    else if (e > -FLINT_BITS)
    {
        s = -e;
        *x2 = d1 >> s;
        *x1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
        *x0 = d0 << (FLINT_BITS - s);
    }
    else if (e == -FLINT_BITS)
    {
        *x1 = d1; *x0 = d0;
    }
    else if (e > -2 * FLINT_BITS)
    {
        s = -e - FLINT_BITS;
        *x1 = d1 >> s;
        *x0 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
    }
    else
        *x0 = d1;
}

/* The reduction of exp: for X = (xi, x2, x1, x0) >= 0 with three fraction
   limbs and xi < 2^62, returns q and sets t = X - q L in [0, L) (three
   fraction limbs), L the floor of log 2 at three limbs: q from
   _mp_real_exp_quotient (exact or one low, corrected once); t is within
   (q + 1) 2^-192 of X - q log 2. */
FLINT_FORCE_INLINE ulong
_arf_small_exp_reduce_3(ulong * t2, ulong * t1, ulong * t0,
    ulong xi, ulong x2, ulong x1, ulong x0)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong q, p3, p2, p1, p0, h, l, u3, u2, u1, u0;

    q = _mp_real_exp_quotient(xi, x2);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(u3, u2, u1, u0, xi, x2, x1, x0, p3, p2, p1, p0);

    if (u3 != 0 || u2 > L[2] || (u2 == L[2] && (u1 > L[1] || (u1 == L[1] && u0 >= L[0]))))
    {
        sub_ddddmmmmssss(u3, u2, u1, u0, u3, u2, u1, u0, UWORD(0), L[2], L[1], L[0]);
        q++;
    }

    *t2 = u2; *t1 = u1; *t0 = u0;
    return q;
}

/* the same with four fraction limbs and L at four limbs (within (q + 1)
   2^-256) */
FLINT_FORCE_INLINE ulong
_arf_small_exp_reduce_4(ulong * t3, ulong * t2, ulong * t1, ulong * t0,
    ulong xi, ulong x3, ulong x2, ulong x1, ulong x0)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong q, p4, p3, p2, p1, p0, h, l, u4, u3, u2, u1, u0;

    q = _mp_real_exp_quotient(xi, x3);
    umul_ppmm(p1, p0, q, L[0]);
    umul_ppmm(h, l, q, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, q, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    sub_dddddmmmmmsssss(u4, u3, u2, u1, u0, xi, x3, x2, x1, x0, p4, p3, p2, p1, p0);

    if (u4 != 0 || u3 > L[3] || (u3 == L[3] && (u2 > L[2] || (u2 == L[2]
            && (u1 > L[1] || (u1 == L[1] && u0 >= L[0]))))))
    {
        sub_dddddmmmmmsssss(u4, u3, u2, u1, u0, u4, u3, u2, u1, u0, UWORD(0), L[3], L[2], L[1], L[0]);
        q++;
    }

    *t3 = u3; *t2 = u2; *t1 = u1; *t0 = u0;
    return q;
}

/* exp(x) at two resp. three limbs from the reduction above and the
   table-driven kernel (mp_real/exp_fp.c explains the error bounds) */
#define ARF_SMALL_EXP_ERR_2 (MP_REAL_SMALL_EXP_ERR_2 + 3)
#define ARF_SMALL_EXP_ERR_3 (MP_REAL_SMALL_EXP_ERR_3 + 3)

/* exp((-1)^sgnbit X) = Y 2^(k - 2 FLINT_BITS), Y = (y2, y1, y0) in
   [B^2, 2 B^2), within ARF_SMALL_EXP_ERR_2 units, for the fixed-point
   X = (xi, x2, x1, x0) (three fraction limbs, xi < 2^62) */
FLINT_FORCE_INLINE slong
_arf_small_exp_fixed_1(ulong * y2, ulong * y1, ulong * y0, ulong xi, ulong x2, ulong x1, ulong x0, int sgnbit)
{
    ulong q, t2, t1, t0;
    slong k;

    q = _arf_small_exp_reduce_3(&t2, &t1, &t0, xi, x2, x1, x0);

    if (sgnbit)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
        sub_dddmmmsss(t2, t1, t0, L[2], L[1], L[0], t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_2(y2, y1, y0, t2, t1);
    return k;
}

/* the same at three limbs: Y = (y3, y2, y1, y0), exp = Y 2^(k - 3 FLINT_BITS),
   for X = (xi, x3, x2, x1, x0) with four fraction limbs */
FLINT_FORCE_INLINE slong
_arf_small_exp_fixed_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong xi, ulong x3, ulong x2, ulong x1, ulong x0, int sgnbit)
{
    ulong q, t3, t2, t1, t0;
    slong k;

    q = _arf_small_exp_reduce_4(&t3, &t2, &t1, &t0, xi, x3, x2, x1, x0);

    if (sgnbit)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
        sub_ddddmmmmssss(t3, t2, t1, t0, L[3], L[2], L[1], L[0], t3, t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_3(y3, y2, y1, y0, t3, t2, t1);
    return k;
}

/* exp(x) = Y 2^(k - 2 FLINT_BITS), Y = (y2, y1, y0) in [B^2, 2 B^2), within
   ARF_SMALL_EXP_ERR_2 units, for the one-limb mantissa d with
   -FLINT_BITS <= e <= MP_REAL_EXP_FP_MAX_EXP */
FLINT_FORCE_INLINE slong
_arf_small_exp_1(ulong * y2, ulong * y1, ulong * y0, ulong d, slong e, int sgnbit)
{
    ulong xi, x1, x0;

    _arf_small_fixed_1(&xi, &x1, &x0, d, e);
    return _arf_small_exp_fixed_1(y2, y1, y0, xi, x1, x0, 0, sgnbit);
}

/* the same at three limbs: Y = (y3, y2, y1, y0), exp(x) = Y 2^(k - 3 FLINT_BITS),
   for -2 FLINT_BITS <= e <= MP_REAL_EXP_FP_MAX_EXP */
FLINT_FORCE_INLINE slong
_arf_small_exp_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0, ulong d1, ulong d0, slong e, int sgnbit)
{
    ulong xi, x2, x1, x0;

    _arf_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);
    return _arf_small_exp_fixed_2(y3, y2, y1, y0, xi, x2, x1, x0, 0, sgnbit);
}

/* |log u| for u = M 2^E, M = (m1, m0) a normalized two-limb fraction and
   |E - 1| < 2^62: c L +- F with c = E - 1, F = log1p(2M - 1) by the
   table-driven kernel at two limbs (2M - 1 exact) and L the floor of
   log 2 at three fraction limbs (c times its error below a quarter ulp of
   B^-2), as (r2, r1, r0) with the integral limb r2, the lowest limb of
   the sum dropped: within ARF_SMALL_LOG_ERR_2 units of B^-2. Returns
   1 if log u < 0 (c < 0). */
#define ARF_SMALL_LOG_ERR_2 (MP_REAL_SMALL_LOG_ERR_2 + 2)
#define ARF_SMALL_LOG_ERR_3 (MP_REAL_SMALL_LOG_ERR_3 + 2)

FLINT_FORCE_INLINE int
_arf_small_log_2(ulong * r2, ulong * r1, ulong * r0, ulong m1, ulong m0, slong E)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
    ulong f1, f0, cc, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_2(&f1, &f0, (m1 << 1) | (m0 >> (FLINT_BITS - 1)), m0 << 1);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    (void) p0;

    if (c >= 0)
        add_sssaaaaaa(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);
    else
        sub_dddmmmsss(p3, p2, p1, p3, p2, p1, UWORD(0), f1, f0);

    *r2 = p3; *r1 = p2; *r0 = p1;
    return c < 0;
}

/* the same at three limbs: M = (m2, m1, m0), the result (r3, r2, r1, r0)
   with three fraction limbs, within ARF_SMALL_LOG_ERR_3 units of B^-3 */
FLINT_FORCE_INLINE int
_arf_small_log_3(ulong * r3, ulong * r2, ulong * r1, ulong * r0,
    ulong m2, ulong m1, ulong m0, slong E)
{
    nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
    ulong f2, f1, f0, cc, p4, p3, p2, p1, p0, h, l;
    slong c = E - 1;

    _mp_real_small_log1p_3(&f2, &f1, &f0, (m2 << 1) | (m1 >> (FLINT_BITS - 1)),
        (m1 << 1) | (m0 >> (FLINT_BITS - 1)), m0 << 1);

    cc = (c < 0) ? -(ulong) c : (ulong) c;
    umul_ppmm(p1, p0, cc, L[0]);
    umul_ppmm(h, l, cc, L[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, cc, L[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, cc, L[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    (void) p0;

    if (c >= 0)
        add_ssssaaaaaaaa(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);
    else
        sub_ddddmmmmssss(p4, p3, p2, p1, p4, p3, p2, p1, UWORD(0), f2, f1, f0);

    *r3 = p4; *r2 = p3; *r1 = p2; *r0 = p1;
    return c < 0;
}

/* whether the register path of log applies to the mantissa (d1 the top
   limb) with exponent E: away from 1 by at least 2^-32 resp. 2^-24 (the
   cancellation would cost too many of the 128 resp. 192 bits) */
#define ARF_SMALL_LOG_OK_1(d1, E) \
    (((E) != 0 && (E) != 1) || ((E) == 1 && ((d1) << 1) >= (UWORD(1) << 32)) \
        || ((E) == 0 && (d1) <= UWORD_MAX - (UWORD(1) << 32)))
#define ARF_SMALL_LOG_OK_2(d1, E) \
    (((E) != 0 && (E) != 1) || ((E) == 1 && ((d1) << 1) >= (UWORD(1) << 25)) \
        || ((E) == 0 && (d1) <= UWORD_MAX - (UWORD(1) << 24)))

/* log1p(x) at two resp. three limbs: (r2, r1, r0) resp. (r3, r2, r1, r0)
   with two resp. three fraction limbs within *err ulps, *neg the sign;
   for 2^-41 <= x < 1 the kernel on x directly, for 1 <= x < 2^63 the
   logarithm of u = 1 + x (exact on two resp. three limbs, no
   cancellation), for -1 < x <= -2^-32 resp. -2^-40 the logarithm of
   u = 1 - |x|, also exact.  Returns 0, or ARF_ELEM_UNCERTAIN for
   other x. */
FLINT_FORCE_INLINE int
_arf_small_log1p_1(ulong * r2, ulong * r1, ulong * r0, ulong * err, int * neg, ulong d, slong e, int sgnbit)
{
    ulong m1, m0, xi, xf;
    slong E;
    unsigned int c;

    if (!sgnbit)
    {
        if (e <= 0)
        {
            if (e < -40)
                return ARF_ELEM_UNCERTAIN;
            _mp_real_small_log1p_2(r1, r0, (e == 0) ? d : (d >> (-e)),
                (e == 0) ? 0 : (d << (FLINT_BITS + e)));
            *r2 = 0;
            *err = MP_REAL_SMALL_LOG_ERR_2;
            *neg = 0;
            return 0;
        }

        if (e >= FLINT_BITS)
            return ARF_ELEM_UNCERTAIN;

        /* u = 1 + x = (xi, xf) with the binary point between the limbs */
        xi = (d >> (FLINT_BITS - e)) + 1;
        xf = d << e;
        c = flint_clz(xi);
        m1 = (c == 0) ? xi : (xi << c) | (xf >> (FLINT_BITS - c));
        m0 = xf << c;
        E = FLINT_BITS - c;
    }
    else
    {
        if (e > 0 || e < -31)
            return ARF_ELEM_UNCERTAIN;

        /* u = 1 - |x|, the top limb nonzero (1 - |x| >= 2^-64 on the
           one-limb mantissa) */
        m1 = (e == 0) ? d : (d >> (-e));
        m0 = (e == 0) ? 0 : (d << (FLINT_BITS + e));
        sub_ddmmss(m1, m0, UWORD(0), UWORD(0), m1, m0);
        c = flint_clz(m1);
        if (c != 0)
        {
            m1 = (m1 << c) | (m0 >> (FLINT_BITS - c));
            m0 = m0 << c;
        }
        E = -(slong) c;
    }

    *neg = _arf_small_log_2(r2, r1, r0, m1, m0, E);
    *err = ARF_SMALL_LOG_ERR_2;
    return 0;
}

FLINT_FORCE_INLINE int
_arf_small_log1p_2(ulong * r3, ulong * r2, ulong * r1, ulong * r0, ulong * err, int * neg,
    ulong d1, ulong d0, slong e, int sgnbit)
{
    ulong m2, m1, m0, xi;
    slong E;
    unsigned int c, s;

    if (!sgnbit)
    {
        if (e <= 0)
        {
            if (e < -40)
                return ARF_ELEM_UNCERTAIN;
            if (e == 0)
                _mp_real_small_log1p_3(r2, r1, r0, d1, d0, 0);
            else
            {
                s = -e;
                _mp_real_small_log1p_3(r2, r1, r0, d1 >> s,
                    (d1 << (FLINT_BITS - s)) | (d0 >> s), d0 << (FLINT_BITS - s));
            }
            *r3 = 0;
            *err = MP_REAL_SMALL_LOG_ERR_3;
            *neg = 0;
            return 0;
        }

        if (e >= FLINT_BITS)
            return ARF_ELEM_UNCERTAIN;

        /* u = 1 + x = (xi, m1, m0) */
        xi = (d1 >> (FLINT_BITS - e)) + 1;
        m1 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        m0 = d0 << e;
        c = flint_clz(xi);
        if (c == 0)
        {
            m2 = xi;
        }
        else
        {
            m2 = (xi << c) | (m1 >> (FLINT_BITS - c));
            m1 = (m1 << c) | (m0 >> (FLINT_BITS - c));
            m0 = m0 << c;
        }
        E = FLINT_BITS - c;
    }
    else
    {
        if (e > 0 || e < -39)
            return ARF_ELEM_UNCERTAIN;

        if (e == 0)
        {
            m2 = d1; m1 = d0; m0 = 0;
        }
        else
        {
            s = -e;
            m2 = d1 >> s;
            m1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
            m0 = d0 << (FLINT_BITS - s);
        }
        sub_dddmmmsss(m2, m1, m0, UWORD(0), UWORD(0), UWORD(0), m2, m1, m0);
        E = 0;
        if (m2 == 0)
        {
            /* |x| within 2^-64 of 1 (e = 0): then m1 != 0 */
            m2 = m1; m1 = m0; m0 = 0;
            E = -FLINT_BITS;
        }
        c = flint_clz(m2);
        if (c != 0)
        {
            m2 = (m2 << c) | (m1 >> (FLINT_BITS - c));
            m1 = (m1 << c) | (m0 >> (FLINT_BITS - c));
            m0 = m0 << c;
        }
        E -= (slong) c;
    }

    *neg = _arf_small_log_3(r3, r2, r1, r0, m2, m1, m0, E);
    *err = ARF_SMALL_LOG_ERR_3;
    return 0;
}

/*
    Register reductions for the trigonometric functions at one and two
    limbs (64-bit only). The reductions give |x| = a pi/2 + s v resp.
    pi |x| = a pi/2 + s v with v in [0, pi/4] as a fixed-point number for
    the table-driven sine and cosine kernel of mp_real at two resp. three
    limbs, within ev ulps of the kernel:

    * x itself for 2^-33 <= |x| < 2^63 (resp. 2^-41): |x| < 1 goes to the
      kernel unreduced (exactly, ev = 0) except for the tangent of
      |x| >= 1/2, larger |x| is reduced with q from a one-limb 2/pi (one
      low at most, corrected once) and 2q times the floor P of pi/4 at
      three resp. four fraction limbs (an error of 2(q + 1) units of B^-3
      resp. B^-4, below 0.65 ulp of B^-2 resp. B^-3 for q < 2^62.4),
      folded to v in [0, pi/4]; the kernel input drops the last limb
      (ev = 2).

    * pi x for 2^-32 <= |x| < 2^63 (resp. 2^-40): 2|x| = m + f exactly,
      t = f/2 or (1 - f)/2 in [0, 1/4] (exact at three limbs, or within
      half an ulp at two limbs), and v = 3t + (pi - 3) t (the product
      truncated, pi - 3 floored: within 6 ulps of B^-3; ev = 2 resp. 8).
      Half-integers are left to the general code (exact results).

    The reductions return ARF_ELEM_UNCERTAIN when v < 2^-32 resp. 2^-40,
    where sin v (or tan v) would lose too much relative accuracy.
*/

#define ARF_SMALL_TRIG_SIGNS(negs_out, negc_out, a, sneg, xsgn) \
    do { \
        int __negs, __negc; \
        if ((a) & 1) { __negs = ((a) == 3); __negc = (sneg) ^ ((a) == 1); } \
        else { __negs = (sneg) ^ ((a) == 2); __negc = ((a) == 2); } \
        /* sin |x| from sin v (a even) or cos v (a odd), and conversely */ \
        negs_out = ((a) & 1) ? __negc : __negs; \
        negc_out = ((a) & 1) ? __negs : __negc; \
    } while (0)

/* reduction of x: v = (v1, v0) at two fraction limbs (ev = 2 for |x| >= 1) */
FLINT_FORCE_INLINE int
_arf_small_trig_reduce_1(int * ra, int * rsneg, ulong * v1, ulong * v0, ulong d, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
    ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
    ulong xi, x1, q, q2, h, l, h2, l2, p3, p2, p1, p0, t3, t2, t1, t0, u3, u2, u1, u0;
    int sneg = 0;

    /* |x| = (xi, x1) with the binary point after xi; 0 <= e < FLINT_BITS */
    xi = (e == 0) ? 0 : (d >> (FLINT_BITS - e));
    x1 = (e == 0) ? d : (d << e);

    /* q = floor(|x| 2/pi), or one less (the one-limb K errs by
       |x| 2^-64 < 1/2) */
    umul_ppmm(h, l, x1, K);
    umul_ppmm(h2, l2, xi, K);
    add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
    (void) l;

    /* t = |x| - 2q P >= 0 */
    q2 = 2 * q;
    umul_ppmm(p1, p0, q2, P[0]);
    umul_ppmm(h, l, q2, P[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q2, P[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(t3, t2, t1, t0, xi, x1, UWORD(0), UWORD(0), p3, p2, p1, p0);

    /* u = 2P (pi/2) */
    u3 = P[2] >> (FLINT_BITS - 1);
    u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
    u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    u0 = P[0] << 1;

    if (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2 && (t1 > u1 || (t1 == u1 && t0 >= u0))))))
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, t3, t2, t1, t0, u3, u2, u1, u0);
        q++;
    }

    /* fold t > pi/4 */
    if (t3 != 0 || t2 > P[2] || (t2 == P[2] && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))
    {
        sub_ddddmmmmssss(t3, t2, t1, t0, u3, u2, u1, u0, t3, t2, t1, t0);
        q++;
        sneg = 1;
    }

    if (t2 < (UWORD(1) << 32))
        return ARF_ELEM_UNCERTAIN;

    *ra = q & 3;
    *rsneg = sneg;
    *v1 = t2;
    *v0 = t1;
    return 0;
}

FLINT_FORCE_INLINE int
_arf_small_trig_reduce_2(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong d1, ulong d0, slong e)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
    ulong K = _mp_real_const_ptr(MP_REAL_CONST_ID_2_DIV_PI, 1)[0];
    ulong xi, x2, x1, q, q2, h, l, h2, l2, p4, p3, p2, p1, p0;
    ulong t4, t3, t2, t1, t0, u4, u3, u2, u1, u0;
    int sneg = 0;

    /* |x| = (xi, x2, x1) with the binary point after xi; 0 <= e < FLINT_BITS */
    if (e == 0)
    {
        xi = 0; x2 = d1; x1 = d0;
    }
    else
    {
        xi = d1 >> (FLINT_BITS - e);
        x2 = (d1 << e) | (d0 >> (FLINT_BITS - e));
        x1 = d0 << e;
    }

    umul_ppmm(h, l, x2, K);
    umul_ppmm(h2, l2, xi, K);
    add_ssaaaa(q, l2, h2, l2, UWORD(0), h);
    (void) l;

    q2 = 2 * q;
    umul_ppmm(p1, p0, q2, P[0]);
    umul_ppmm(h, l, q2, P[1]);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, q2, P[2]);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    umul_ppmm(h, l, q2, P[3]);
    add_ssaaaa(p4, p3, h, l, UWORD(0), p3);
    sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, xi, x2, x1, UWORD(0), UWORD(0), p4, p3, p2, p1, p0);

    u4 = P[3] >> (FLINT_BITS - 1);
    u3 = (P[3] << 1) | (P[2] >> (FLINT_BITS - 1));
    u2 = (P[2] << 1) | (P[1] >> (FLINT_BITS - 1));
    u1 = (P[1] << 1) | (P[0] >> (FLINT_BITS - 1));
    u0 = P[0] << 1;

    if (t4 > u4 || (t4 == u4 && (t3 > u3 || (t3 == u3 && (t2 > u2 || (t2 == u2
            && (t1 > u1 || (t1 == u1 && t0 >= u0))))))))
    {
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, t4, t3, t2, t1, t0, u4, u3, u2, u1, u0);
        q++;
    }

    if (t4 != 0 || t3 > P[3] || (t3 == P[3] && (t2 > P[2] || (t2 == P[2]
            && (t1 > P[1] || (t1 == P[1] && t0 > P[0]))))))
    {
        sub_dddddmmmmmsssss(t4, t3, t2, t1, t0, u4, u3, u2, u1, u0, t4, t3, t2, t1, t0);
        q++;
        sneg = 1;
    }

    if (t3 < (UWORD(1) << 24))
        return ARF_ELEM_UNCERTAIN;

    *ra = q & 3;
    *rsneg = sneg;
    *v2 = t3;
    *v1 = t2;
    *v0 = t1;
    return 0;
}

/* v = pi t = 3t + (pi - 3) t for t = (t2, t1, t0) <= 1/4 */
FLINT_FORCE_INLINE void
_arf_small_trig_pi_times(ulong * v2, ulong * v1, ulong * v0, ulong t2, ulong t1, ulong t0)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 4);
    ulong c2, c1, c0, s2, s1, s0;

    /* pi - 3 at three fraction limbs, floored: the fraction of 4P */
    c2 = (P[3] << 2) | (P[2] >> (FLINT_BITS - 2));
    c1 = (P[2] << 2) | (P[1] >> (FLINT_BITS - 2));
    c0 = (P[1] << 2) | (P[0] >> (FLINT_BITS - 2));

    _mp_real_small_mulhi_3x3(&s2, &s1, &s0, c2, c1, c0, t2, t1, t0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, t2, t1, t0);
    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, (t2 << 1) | (t1 >> (FLINT_BITS - 1)),
        (t1 << 1) | (t0 >> (FLINT_BITS - 1)), t0 << 1);
    *v2 = s2; *v1 = s1; *v0 = s0;
}

/* 2|x| = m + (f2, f1, f0) to the reduction pi |x| = a pi/2 + s v: t at
   three limbs (exact if f0's lowest bit is zero) */
FLINT_FORCE_INLINE int
_arf_small_trig_pi_fold(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong m, ulong f2, ulong f1, ulong f0)
{
    if (f2 == 0 && f1 == 0 && f0 == 0)
        return ARF_ELEM_UNCERTAIN;

    if (f2 < (UWORD(1) << (FLINT_BITS - 1)) || (f2 == (UWORD(1) << (FLINT_BITS - 1)) && f1 == 0 && f0 == 0))
    {
        *ra = m & 3;
        *rsneg = 0;
    }
    else
    {
        *ra = (m + 1) & 3;
        *rsneg = 1;
        sub_dddmmmsss(f2, f1, f0, UWORD(0), UWORD(0), UWORD(0), f2, f1, f0);
    }

    /* t = f / 2 */
    _arf_small_trig_pi_times(v2, v1, v0, f2 >> 1, (f2 << (FLINT_BITS - 1)) | (f1 >> 1),
        (f1 << (FLINT_BITS - 1)) | (f0 >> 1));
    return 0;
}

/* reduction of pi x, one limb, -31 <= e < FLINT_BITS */
FLINT_FORCE_INLINE int
_arf_small_trig_pi_reduce_1(int * ra, int * rsneg, ulong * v1, ulong * v0, ulong d, slong e)
{
    ulong m, f1, f0, v2, w1, w0;

    if (e >= 1)
    {
        if (e == FLINT_BITS - 1)
            return ARF_ELEM_UNCERTAIN;
        m = d >> (FLINT_BITS - 1 - e);
        f1 = d << (e + 1);
        f0 = 0;
    }
    else if (e == 0)
    {
        m = 1;
        f1 = d << 1;
        f0 = 0;
    }
    else
    {
        unsigned int s = -e - 1;
        m = 0;
        f1 = d >> s;
        f0 = (s == 0) ? 0 : (d << (FLINT_BITS - s));
    }

    if (_arf_small_trig_pi_fold(ra, rsneg, &v2, &w1, &w0, m, f1, f0, 0) != 0)
        return ARF_ELEM_UNCERTAIN;

    if (v2 < (UWORD(1) << 32))
        return ARF_ELEM_UNCERTAIN;

    *v1 = v2;
    *v0 = w1;
    return 0;
}

/* reduction of pi x, two limbs, -39 <= e < FLINT_BITS */
FLINT_FORCE_INLINE int
_arf_small_trig_pi_reduce_2(int * ra, int * rsneg, ulong * v2, ulong * v1, ulong * v0,
    ulong d1, ulong d0, slong e)
{
    ulong m, f2, f1, f0;

    if (e >= 1)
    {
        if (e == FLINT_BITS - 1)
        {
            m = d1;
            f2 = d0;
            f1 = 0;
        }
        else
        {
            m = d1 >> (FLINT_BITS - 1 - e);
            f2 = (d1 << (e + 1)) | (d0 >> (FLINT_BITS - 1 - e));
            f1 = d0 << (e + 1);
        }
        f0 = 0;
    }
    else if (e == 0)
    {
        m = 1;
        f2 = (d1 << 1) | (d0 >> (FLINT_BITS - 1));
        f1 = d0 << 1;
        f0 = 0;
    }
    else
    {
        unsigned int s = -e - 1;
        m = 0;
        if (s == 0)
        {
            f2 = d1; f1 = d0; f0 = 0;
        }
        else
        {
            f2 = d1 >> s;
            f1 = (d1 << (FLINT_BITS - s)) | (d0 >> s);
            f0 = d0 << (FLINT_BITS - s);
        }
    }

    if (_arf_small_trig_pi_fold(ra, rsneg, v2, v1, v0, m, f2, f1, f0) != 0)
        return ARF_ELEM_UNCERTAIN;

    if (*v2 < (UWORD(1) << 24))
        return ARF_ELEM_UNCERTAIN;

    return 0;
}

/* (y2, y1, y0) >> s for 0 <= s, truncated (zero for s >= 3 FLINT_BITS) */
FLINT_FORCE_INLINE void
_arf_small_rshift_3(ulong * r2, ulong * r1, ulong * r0, ulong y2, ulong y1, ulong y0, ulong s)
{
    unsigned int b;

    if (s >= 2 * FLINT_BITS)
    {
        *r2 = 0; *r1 = 0;
        *r0 = (s >= 3 * FLINT_BITS) ? 0 : y2 >> (s - 2 * FLINT_BITS);
        return;
    }
    if (s >= FLINT_BITS)
    {
        y0 = y1; y1 = y2; y2 = 0;
        s -= FLINT_BITS;
    }
    b = s;
    if (b != 0)
    {
        y0 = (y0 >> b) | (y1 << (FLINT_BITS - b));
        y1 = (y1 >> b) | (y2 << (FLINT_BITS - b));
        y2 = y2 >> b;
    }
    *r2 = y2; *r1 = y1; *r0 = y0;
}

/* (y3, y2, y1, y0) >> s, truncated */
FLINT_FORCE_INLINE void
_arf_small_rshift_4(ulong * r3, ulong * r2, ulong * r1, ulong * r0,
    ulong y3, ulong y2, ulong y1, ulong y0, ulong s)
{
    unsigned int b;

    if (s >= 4 * FLINT_BITS)
    {
        *r3 = *r2 = *r1 = *r0 = 0;
        return;
    }
    while (s >= FLINT_BITS)
    {
        y0 = y1; y1 = y2; y2 = y3; y3 = 0;
        s -= FLINT_BITS;
    }
    b = s;
    if (b != 0)
    {
        y0 = (y0 >> b) | (y1 << (FLINT_BITS - b));
        y1 = (y1 >> b) | (y2 << (FLINT_BITS - b));
        y2 = (y2 >> b) | (y3 << (FLINT_BITS - b));
        y3 = y3 >> b;
    }
    *r3 = y3; *r2 = y2; *r1 = y1; *r0 = y0;
}
/* expm1(x) at two resp. three limbs for -40 <= e <= MP_REAL_EXP_FP_MAX_EXP
   (the cancellation below 2^-41 would cost too many of the 128 resp.
   192 bits): expm1(x) = (-1)^(*sgn) Y 2^(k - 2 FLINT_BITS) with
   Y = (y2, y1, y0) within ARF_SMALL_EXP_ERR_2 + 1 units, from
   exp(x) = Y' 2^(k' - 2 FLINT_BITS) by subtracting 1 in the same fixed
   point (2^k' Y' - 1 = 2^k' (Y' - 2^-k') for x > 0, 1 - 2^-(q+1) Y' for
   x < 0) */
FLINT_FORCE_INLINE slong
_arf_small_expm1_1(ulong * y2, ulong * y1, ulong * y0, int * sgn, ulong d, slong e, int sgnbit)
{
    slong k = _arf_small_exp_1(y2, y1, y0, d, e, sgnbit);

    if (!sgnbit)
    {
        /* Y - 2^-k, k >= 0: the bit 2 FLINT_BITS - k of Y (or below an
           ulp) */
        if (k <= 2 * FLINT_BITS)
        {
            slong b = 2 * FLINT_BITS - k;
            ulong b2 = 0, b1 = 0, b0 = 0;
            if (b >= 2 * FLINT_BITS)
                b2 = UWORD(1);
            else if (b >= FLINT_BITS)
                b1 = UWORD(1) << (b - FLINT_BITS);
            else
                b0 = UWORD(1) << b;
            sub_dddmmmsss(*y2, *y1, *y0, *y2, *y1, *y0, b2, b1, b0);
        }
        *sgn = 0;
        return k;
    }
    else
    {
        /* 1 - Y 2^(k - 2 FLINT_BITS), k <= -1 (one more unit for the
           shifted-out bits) */
        _arf_small_rshift_3(y2, y1, y0, *y2, *y1, *y0, -k);
        sub_dddmmmsss(*y2, *y1, *y0, UWORD(1), UWORD(0), UWORD(0), *y2, *y1, *y0);
        *sgn = 1;
        return 0;
    }
}

FLINT_FORCE_INLINE slong
_arf_small_expm1_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0, int * sgn, ulong d1, ulong d0, slong e, int sgnbit)
{
    slong k = _arf_small_exp_2(y3, y2, y1, y0, d1, d0, e, sgnbit);

    if (!sgnbit)
    {
        if (k <= 3 * FLINT_BITS)
        {
            slong b = 3 * FLINT_BITS - k;
            ulong b3 = 0, b2 = 0, b1 = 0, b0 = 0;
            if (b >= 3 * FLINT_BITS)
                b3 = UWORD(1);
            else if (b >= 2 * FLINT_BITS)
                b2 = UWORD(1) << (b - 2 * FLINT_BITS);
            else if (b >= FLINT_BITS)
                b1 = UWORD(1) << (b - FLINT_BITS);
            else
                b0 = UWORD(1) << b;
            sub_ddddmmmmssss(*y3, *y2, *y1, *y0, *y3, *y2, *y1, *y0, b3, b2, b1, b0);
        }
        *sgn = 0;
        return k;
    }
    else
    {
        _arf_small_rshift_4(y3, y2, y1, y0, *y3, *y2, *y1, *y0, -k);
        sub_ddddmmmmssss(*y3, *y2, *y1, *y0, UWORD(1), UWORD(0), UWORD(0), UWORD(0), *y3, *y2, *y1, *y0);
        *sgn = 1;
        return 0;
    }
}

/* Arctangents at two resp. three limbs (64-bit), with the table-driven
   kernels of mp_real at two resp. three fraction limbs:

   atan |x| for 2^-33 <= |x| < 1 (resp. 2^-41) by the kernel on |x|
   exactly; for |x| >= 1 with e <= 2 FLINT_BITS (resp. 3 FLINT_BITS), as
   pi/2 - atan(Y / X) with Y = 2^-e and X = M exactly (one division in
   the ratio kernel).

   atan2 as c pi/4 + sigma atan(Y / X) with Y = Mnum 2^-d (exact at
   two resp. three limbs for d < FLINT_BITS) and X = Mden, for d <= 32
   (resp. 40) when c = 0, where atan(t) >= 2^-34 must keep its relative
   accuracy, and d < FLINT_BITS otherwise (the result is then at least
   pi/4).

   c pi/4 comes from the floor of pi/4 at two resp. three limbs (c ulps
   below).  The results are frames (t2, t1, t0) resp. (t3, t2, t1, t0)
   with an integral top limb, within *err units, for _arf_small_set_1
   resp. _2 with k = 0; the helpers return ARF_ELEM_UNCERTAIN for the
   other arguments. */

/* the frame c pi/4 + sigma a for a = (a1, a0) within ea units, c in
   {0, 2, 4}, not cancelling */
FLINT_FORCE_INLINE void
_arf_small_atan_finish_1(ulong * t2, ulong * t1, ulong * t0, ulong * err,
    ulong a1, ulong a0, ulong ea, int c, int sigma)
{
    nn_srcptr P;
    unsigned int s;

    if (c == 0)
    {
        *t2 = 0; *t1 = a1; *t0 = a0;
        *err = ea;
        return;
    }

    P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 2);
    s = (c == 2) ? 1 : 2;
    *t2 = P[1] >> (FLINT_BITS - s);
    *t1 = (P[1] << s) | (P[0] >> (FLINT_BITS - s));
    *t0 = P[0] << s;
    if (sigma > 0)
        add_sssaaaaaa(*t2, *t1, *t0, *t2, *t1, *t0, UWORD(0), a1, a0);
    else
        sub_dddmmmsss(*t2, *t1, *t0, *t2, *t1, *t0, UWORD(0), a1, a0);
    *err = ea + c;
}

FLINT_FORCE_INLINE void
_arf_small_atan_finish_2(ulong * t3, ulong * t2, ulong * t1, ulong * t0, ulong * err,
    ulong a2, ulong a1, ulong a0, ulong ea, int c, int sigma)
{
    nn_srcptr P;
    unsigned int s;

    if (c == 0)
    {
        *t3 = 0; *t2 = a2; *t1 = a1; *t0 = a0;
        *err = ea;
        return;
    }

    P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, 3);
    s = (c == 2) ? 1 : 2;
    *t3 = P[2] >> (FLINT_BITS - s);
    *t2 = (P[2] << s) | (P[1] >> (FLINT_BITS - s));
    *t1 = (P[1] << s) | (P[0] >> (FLINT_BITS - s));
    *t0 = P[0] << s;
    if (sigma > 0)
        add_ssssaaaaaaaa(*t3, *t2, *t1, *t0, *t3, *t2, *t1, *t0, UWORD(0), a2, a1, a0);
    else
        sub_ddddmmmmssss(*t3, *t2, *t1, *t0, *t3, *t2, *t1, *t0, UWORD(0), a2, a1, a0);
    *err = ea + c;
}

/* atan |x| for |x| = d 2^(e - FLINT_BITS) */
FLINT_FORCE_INLINE int
_arf_small_atan_1(ulong * t2, ulong * t1, ulong * t0, ulong * err, ulong d, slong e)
{
    ulong a1, a0;

    if (e <= 0)
    {
        if (e < -32)
            return ARF_ELEM_UNCERTAIN;
        _mp_real_small_atan_2(&a1, &a0, (e == 0) ? d : (d >> (-e)),
            (e == 0) ? 0 : (d << (FLINT_BITS + e)));
        _arf_small_atan_finish_1(t2, t1, t0, err, a1, a0, MP_REAL_SMALL_ATAN_ERR_2, 0, 1);
        return 0;
    }

    if (e > 2 * FLINT_BITS)
        return ARF_ELEM_UNCERTAIN;

    /* atan(1/|x|) = atan(2^-e / M) */
    if (e <= FLINT_BITS)
        _mp_real_small_atan_ratio_2(&a1, &a0, UWORD(1) << (FLINT_BITS - e), 0, d, 0);
    else
        _mp_real_small_atan_ratio_2(&a1, &a0, 0, UWORD(1) << (2 * FLINT_BITS - e), d, 0);

    _arf_small_atan_finish_1(t2, t1, t0, err, a1, a0, MP_REAL_SMALL_ATAN_ERR_2, 2, -1);
    return 0;
}

/* atan |x| for |x| = (d1, d0) 2^(e - 2 FLINT_BITS) */
FLINT_FORCE_INLINE int
_arf_small_atan_2(ulong * t3, ulong * t2, ulong * t1, ulong * t0, ulong * err, ulong d1, ulong d0, slong e)
{
    ulong a2, a1, a0;

    if (e <= 0)
    {
        ulong v2, v1, v0;

        if (e < -40)
            return ARF_ELEM_UNCERTAIN;
        if (e == 0)
        {
            v2 = d1; v1 = d0; v0 = 0;
        }
        else
        {
            unsigned int sh = -e;
            v2 = d1 >> sh;
            v1 = (d1 << (FLINT_BITS - sh)) | (d0 >> sh);
            v0 = d0 << (FLINT_BITS - sh);
        }
        _mp_real_small_atan_3(&a2, &a1, &a0, v2, v1, v0);
        _arf_small_atan_finish_2(t3, t2, t1, t0, err, a2, a1, a0, MP_REAL_SMALL_ATAN_ERR_3, 0, 1);
        return 0;
    }

    if (e > 3 * FLINT_BITS)
        return ARF_ELEM_UNCERTAIN;

    if (e <= FLINT_BITS)
        _mp_real_small_atan_ratio_3(&a2, &a1, &a0, UWORD(1) << (FLINT_BITS - e), 0, 0, d1, d0, 0);
    else if (e <= 2 * FLINT_BITS)
        _mp_real_small_atan_ratio_3(&a2, &a1, &a0, 0, UWORD(1) << (2 * FLINT_BITS - e), 0, d1, d0, 0);
    else
        _mp_real_small_atan_ratio_3(&a2, &a1, &a0, 0, 0, UWORD(1) << (3 * FLINT_BITS - e), d1, d0, 0);

    _arf_small_atan_finish_2(t3, t2, t1, t0, err, a2, a1, a0, MP_REAL_SMALL_ATAN_ERR_3, 2, -1);
    return 0;
}

/* atan2 with t = (Mnum / Mden) 2^-d for one-limb mantissas */
FLINT_FORCE_INLINE int
_arf_small_atan2_1(ulong * t2, ulong * t1, ulong * t0, ulong * err,
    ulong num, ulong den, slong d, int c, int sigma)
{
    ulong a1, a0, y1, y0;

    if (d >= FLINT_BITS || (c == 0 && d > 32))
        return ARF_ELEM_UNCERTAIN;

    y1 = num >> d;
    y0 = (d == 0) ? 0 : (num << (FLINT_BITS - d));
    _mp_real_small_atan_ratio_2(&a1, &a0, y1, y0, den, 0);
    _arf_small_atan_finish_1(t2, t1, t0, err, a1, a0, MP_REAL_SMALL_ATAN_ERR_2, c, sigma);
    return 0;
}

/* the same for two-limb mantissas */
FLINT_FORCE_INLINE int
_arf_small_atan2_2(ulong * t3, ulong * t2, ulong * t1, ulong * t0, ulong * err,
    ulong num1, ulong num0, ulong den1, ulong den0, slong d, int c, int sigma)
{
    ulong a2, a1, a0, y2, y1, y0;

    if (d >= FLINT_BITS || (c == 0 && d > 40))
        return ARF_ELEM_UNCERTAIN;

    if (d == 0)
    {
        y2 = num1; y1 = num0; y0 = 0;
    }
    else
    {
        y2 = num1 >> d;
        y1 = (num1 << (FLINT_BITS - d)) | (num0 >> d);
        y0 = num0 << (FLINT_BITS - d);
    }
    _mp_real_small_atan_ratio_3(&a2, &a1, &a0, y2, y1, y0, den1, den0, 0);
    _arf_small_atan_finish_2(t3, t2, t1, t0, err, a2, a1, a0, MP_REAL_SMALL_ATAN_ERR_3, c, sigma);
    return 0;
}

/* Hyperbolic functions at one and two limbs (64-bit): with |x| =
   q log 2 + t (the reduction of the exp fast paths), exp(|x|) = 2^q Y
   and exp(-|x|) = 2^-q Z, Y = exp(t) and Z = exp(-t) by the paired
   kernel _mp_real_small_exp_pm_2 resp. _3 (one evaluation sharing the
   polynomial, about 1.5 times one exp kernel), so that

       cosh x, |sinh x| = 2^(q-1) (Y +- Z 2^-2q)

   in fixed point.  Y is within MP_REAL_SMALL_EXP_ERR units (the kernel,
   the reduction error of t and the dropped limb, times exp' < 2), Z <= 1
   within the kernel's ERR plus the error of t (below 1.4 units), and
   shifted down by 2q bits within one more unit (or dropped when below
   one unit): ERR_ND = 2 ERR + 6 units in all.  Without cancellation for
   cosh; sinh keeps 2^-90 resp. 2^-144 relative accuracy for |x| >=
   2^-32 resp. 2^-41 (the callers' ranges), and e <= MP_REAL_EXP_FP_MAX_EXP.

   tanh and coth from the same pair: with N = Y - Z 2^-2q and
   D = Y + Z 2^-2q in fixed point (each within ERR_ND units), |tanh x| =
   N / D < 1 and |coth x| = D / N > 1 by one division in registers.  For
   tanh, N and D are scaled by 2^-sh into [1/2, 1) (D in [1, 3); one more
   unit for each truncation): within 2 ERR_ND + 4 + 1 units, the error of
   D entering at most with N / D < 1 and D >= 1.  For coth, N = Nn 2^-s
   with Nn in [1/2, 1) and the quotient (D / 8) / Nn at most 3/4: within
   ERR_ND (2^(max(s, 0) + 1) + 1) + 4 units of the quotient, the relative
   error of N being at most 2^(s+1) ERR_ND units.  When Z 2^-2q falls
   below a unit, tanh and coth are 1 within one unit (the helpers return
   ARF_SMALL_HYP_ONE; the callers of nearly-exact precision handle those
   arguments before).

   sech and csch (which = 2, 3) are 2^(1-q) / D resp. 2^(1-q) / N, as
   (1/2 - B^-w) / Dn resp. (1/2 - B^-w) / Nn for the normalized Dn =
   D 2^-sh and Nn = N 2^s: within ERR_ND + 4 resp. ERR_ND 2^(max(s, 0) + 1)
   + 4 units of the quotient.

   The results are frames for _arf_small_set_1 resp. _2 with the shift k. */

#define ARF_SMALL_HYP_ONE 1

/* (the codes which, ARF_SMALL_HYP_TANH etc., precede the 64-bit section) */

/* Y and Z 2^-2q at two fraction limbs for the fixed-point X = (xi, x2,
   x1, x0) (three fraction limbs, xi < 2^62), returning q */
FLINT_FORCE_INLINE ulong
_arf_small_exp_pair_fixed_1(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0,
    ulong xi, ulong x2, ulong x1, ulong x0)
{
    ulong q, t2, t1, t0;

    q = _arf_small_exp_reduce_3(&t2, &t1, &t0, xi, x2, x1, x0);
    if (2 * q < 2 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_2(y2, y1, y0, z2, z1, z0, t2, t1);
        _arf_small_rshift_3(z2, z1, z0, *z2, *z1, *z0, 2 * q);
    }
    else
    {
        _mp_real_small_exp_2(y2, y1, y0, t2, t1);
        *z2 = *z1 = *z0 = 0;
    }
    return q;
}

FLINT_FORCE_INLINE ulong
_arf_small_exp_pair_fixed_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong xi, ulong x3, ulong x2, ulong x1, ulong x0)
{
    ulong q, t3, t2, t1, t0;

    q = _arf_small_exp_reduce_4(&t3, &t2, &t1, &t0, xi, x3, x2, x1, x0);
    if (2 * q < 3 * FLINT_BITS + 2)
    {
        _mp_real_small_exp_pm_3(y3, y2, y1, y0, z3, z2, z1, z0, t3, t2, t1);
        _arf_small_rshift_4(z3, z2, z1, z0, *z3, *z2, *z1, *z0, 2 * q);
    }
    else
    {
        _mp_real_small_exp_3(y3, y2, y1, y0, t3, t2, t1);
        *z3 = *z2 = *z1 = *z0 = 0;
    }
    return q;
}

/* the same for |x| = d 2^(e - FLINT_BITS) resp. (d1, d0) 2^(e - 2 FLINT_BITS) */
FLINT_FORCE_INLINE ulong
_arf_small_exp_pair_1(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0,
    ulong d, slong e)
{
    ulong xi, x1, x0;

    _arf_small_fixed_1(&xi, &x1, &x0, d, e);
    return _arf_small_exp_pair_fixed_1(y2, y1, y0, z2, z1, z0, xi, x1, x0, 0);
}

FLINT_FORCE_INLINE ulong
_arf_small_exp_pair_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong d1, ulong d0, slong e)
{
    ulong xi, x2, x1, x0;

    _arf_small_fixed_2(&xi, &x2, &x1, &x0, d1, d0, e);
    return _arf_small_exp_pair_fixed_2(y3, y2, y1, y0, z3, z2, z1, z0, xi, x2, x1, x0, 0);
}

#define ARF_SMALL_HYP_ERR_2 (2 * MP_REAL_SMALL_EXP_ERR_2 + 6)
#define ARF_SMALL_HYP_ERR_3 (2 * MP_REAL_SMALL_EXP_ERR_3 + 6)

/* |sinh x| = S 2^(k - 2 FLINT_BITS), cosh x = C 2^(k - 2 FLINT_BITS) for
   the frames S = (s2, s1, s0), C = (c2, c1, c0) within ARF_SMALL_HYP_ERR_2
   units (either output NULL); -31 <= e <= MP_REAL_EXP_FP_MAX_EXP */
FLINT_FORCE_INLINE void
_arf_small_sinh_cosh_1(ulong * s2, ulong * s1, ulong * s0, ulong * c2, ulong * c1, ulong * c0,
    slong * k, ulong d, slong e)
{
    ulong q, y2, y1, y0, z2, z1, z0;

    q = _arf_small_exp_pair_1(&y2, &y1, &y0, &z2, &z1, &z0, d, e);
    *k = (slong) q - 1;
    if (s2 != NULL)
        sub_dddmmmsss(*s2, *s1, *s0, y2, y1, y0, z2, z1, z0);
    if (c2 != NULL)
        add_sssaaaaaa(*c2, *c1, *c0, y2, y1, y0, z2, z1, z0);
}

/* the same at three fraction limbs, -40 <= e */
FLINT_FORCE_INLINE void
_arf_small_sinh_cosh_2(ulong * s3, ulong * s2, ulong * s1, ulong * s0,
    ulong * c3, ulong * c2, ulong * c1, ulong * c0, slong * k, ulong d1, ulong d0, slong e)
{
    ulong q, y3, y2, y1, y0, z3, z2, z1, z0;

    q = _arf_small_exp_pair_2(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, d1, d0, e);
    *k = (slong) q - 1;
    if (s3 != NULL)
        sub_ddddmmmmssss(*s3, *s2, *s1, *s0, y3, y2, y1, y0, z3, z2, z1, z0);
    if (c3 != NULL)
        add_ssssaaaaaaaa(*c3, *c2, *c1, *c0, y3, y2, y1, y0, z3, z2, z1, z0);
}

/* tanh, coth, sech or csch (which) of |x|: the frame (t2, t1, t0)
   2^(k - 2 FLINT_BITS) within *err units; returns ARF_SMALL_HYP_ONE
   (nothing set) when tanh or coth is 1 within a unit */
FLINT_FORCE_INLINE int
_arf_small_tanh_1(ulong * t2, ulong * t1, ulong * t0, slong * k, ulong * err, ulong d, slong e, int which)
{
    ulong y2, y1, y0, z2, z1, z0, n2, n1, n0, d2, d1, d0, q[2], qe;
    ulong e_nd = ARF_SMALL_HYP_ERR_2;
    unsigned int sh;
    slong s;

    qe = _arf_small_exp_pair_1(&y2, &y1, &y0, &z2, &z1, &z0, d, e);

    if (which <= ARF_SMALL_HYP_COTH && z2 == 0 && z1 == 0 && z0 == 0)
        return ARF_SMALL_HYP_ONE;

    sub_dddmmmsss(n2, n1, n0, y2, y1, y0, z2, z1, z0);
    add_sssaaaaaa(d2, d1, d0, y2, y1, y0, z2, z1, z0);

    if (which == ARF_SMALL_HYP_SECH)
    {
        /* sech = 2^(1-q) / D */
        sh = (d2 >= 2) ? 2 : 1;
        _arf_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, sh);
        _mp_real_divq_4_2(q, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX, d1, d0);
        *t2 = 0; *t1 = q[1]; *t0 = q[0];
        *k = 2 - (slong) qe - (slong) sh;
        *err = e_nd + 4;
        return 0;
    }

    if (which == ARF_SMALL_HYP_TANH)
    {
        /* N / D with both scaled into fractions, D normalized */
        sh = (d2 >= 2) ? 2 : 1;
        _arf_small_rshift_3(&n2, &n1, &n0, n2, n1, n0, sh);
        _arf_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, sh);
        _mp_real_divq_4_2(q, n1, n0, UWORD(0), UWORD(0), d1, d0);
        *t2 = 0; *t1 = q[1]; *t0 = q[0];
        *k = 0;
        *err = 2 * e_nd + 5;
        return 0;
    }

    /* N = Nn 2^-s with Nn in [1/2, 1) */
    if (n2 != 0)
    {
        n0 = (n0 >> 1) | (n1 << (FLINT_BITS - 1));
        n1 = (n1 >> 1) | (n2 << (FLINT_BITS - 1));
        s = -1;
    }
    else
    {
        sh = flint_clz(n1);
        if (sh != 0)
        {
            n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
            n0 = n0 << sh;
        }
        s = sh;
    }

    if (which == ARF_SMALL_HYP_CSCH)
    {
        /* csch = 2^(1-q) / N */
        _mp_real_divq_4_2(q, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, UWORD_MAX, n1, n0);
        *t2 = 0; *t1 = q[1]; *t0 = q[0];
        *k = 2 - (slong) qe + s;
        *err = e_nd * (UWORD(2) << FLINT_MAX(s, 0)) + 4;
        return 0;
    }

    /* coth = (D / 8) / Nn */
    _arf_small_rshift_3(&d2, &d1, &d0, d2, d1, d0, 3);
    _mp_real_divq_4_2(q, d1, d0, UWORD(0), UWORD(0), n1, n0);
    *t2 = 0; *t1 = q[1]; *t0 = q[0];
    *k = s + 3;
    *err = e_nd * ((UWORD(2) << FLINT_MAX(s, 0)) + 1) + 4;
    return 0;
}

FLINT_FORCE_INLINE int
_arf_small_tanh_2(ulong * t3, ulong * t2, ulong * t1, ulong * t0, slong * k, ulong * err,
    ulong dd1, ulong dd0, slong e, int which)
{
    ulong y3, y2, y1, y0, z3, z2, z1, z0, n3, n2, n1, n0, d3, d2, d1, d0, w2, w1, w0, qe;
    ulong e_nd = ARF_SMALL_HYP_ERR_3;
    unsigned int sh;
    slong s;

    qe = _arf_small_exp_pair_2(&y3, &y2, &y1, &y0, &z3, &z2, &z1, &z0, dd1, dd0, e);

    if (which <= ARF_SMALL_HYP_COTH && z3 == 0 && z2 == 0 && z1 == 0 && z0 == 0)
        return ARF_SMALL_HYP_ONE;

    sub_ddddmmmmssss(n3, n2, n1, n0, y3, y2, y1, y0, z3, z2, z1, z0);
    add_ssssaaaaaaaa(d3, d2, d1, d0, y3, y2, y1, y0, z3, z2, z1, z0);

    if (which == ARF_SMALL_HYP_SECH)
    {
        sh = (d3 >= 2) ? 2 : 1;
        _arf_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, sh);
        _mp_real_small_divq_3(&w2, &w1, &w0, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, d2, d1, d0);
        *t3 = 0; *t2 = w2; *t1 = w1; *t0 = w0;
        *k = 2 - (slong) qe - (slong) sh;
        *err = e_nd + 3 + MP_REAL_SMALL_DIVQ_3_ERR;
        return 0;
    }

    if (which == ARF_SMALL_HYP_TANH)
    {
        sh = (d3 >= 2) ? 2 : 1;
        _arf_small_rshift_4(&n3, &n2, &n1, &n0, n3, n2, n1, n0, sh);
        _arf_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, sh);
        _mp_real_small_divq_3(&w2, &w1, &w0, n2, n1, n0, d2, d1, d0);
        *t3 = 0; *t2 = w2; *t1 = w1; *t0 = w0;
        *k = 0;
        *err = 2 * e_nd + 4 + MP_REAL_SMALL_DIVQ_3_ERR;
        return 0;
    }

    if (n3 != 0)
    {
        _arf_small_rshift_4(&n3, &n2, &n1, &n0, n3, n2, n1, n0, 1);
        s = -1;
    }
    else
    {
        sh = flint_clz(n2);
        if (sh != 0)
        {
            n2 = (n2 << sh) | (n1 >> (FLINT_BITS - sh));
            n1 = (n1 << sh) | (n0 >> (FLINT_BITS - sh));
            n0 = n0 << sh;
        }
        s = sh;
    }

    if (which == ARF_SMALL_HYP_CSCH)
    {
        _mp_real_small_divq_3(&w2, &w1, &w0, (UWORD(1) << (FLINT_BITS - 1)) - 1, UWORD_MAX, UWORD_MAX, n2, n1, n0);
        *t3 = 0; *t2 = w2; *t1 = w1; *t0 = w0;
        *k = 2 - (slong) qe + s;
        *err = e_nd * (UWORD(2) << FLINT_MAX(s, 0)) + 3 + MP_REAL_SMALL_DIVQ_3_ERR;
        return 0;
    }

    _arf_small_rshift_4(&d3, &d2, &d1, &d0, d3, d2, d1, d0, 3);
    _mp_real_small_divq_3(&w2, &w1, &w0, d2, d1, d0, n2, n1, n0);
    *t3 = 0; *t2 = w2; *t1 = w1; *t0 = w0;
    *k = s + 3;
    *err = e_nd * ((UWORD(2) << FLINT_MAX(s, 0)) + 1) + 3 + MP_REAL_SMALL_DIVQ_3_ERR;
    return 0;
}

/*
    Powers x^y = exp(y log |x|) in registers at one and two limbs (64-bit
    only), for the mantissas (dx, ex), (dy, ey) of |x| and |y| (returning
    the result Y 2^(k - 2 resp. 3 FLINT_BITS) within err units, or
    ARF_ELEM_UNCERTAIN when the path does not apply): log |x| =
    Lambda at three fraction limbs (_arf_small_log_3 on the mantissa,
    within ARF_SMALL_LOG_ERR_3 = 26 units of B^-3, also near |x| = 1
    since only the absolute error of t matters; at one limb and |y| <
    2^40 the two-limb version within 18 units of B^-2 suffices, 2^-82
    relative to the result), t = y Lambda exactly
    from the mantissas and truncated to three resp. four fraction limbs,
    then the reduction and kernel of the exp fast paths. The error of t,
    below (26 |y| + 1) 2^-192, adds twice as many units of the result
    (exp(t +- delta) = exp(t)(1 +- 2 delta)): within 28 units of B^-2 for
    |y| < 2^63 at one limb, and below 2^-146 relative for |y| < 2^40 at
    two limbs. Other y, and |t| >= 2^61 (beyond the exponent range), take
    the general code.
*/

/* (r[0], ..., r[len - 1]) = (a[0], ..., a[len - 1]) >> s for 1 <= s <
   2 FLINT_BITS, truncated */
/* r = a >> s for 0 <= s < 2 FLINT_BITS (len limbs, truncated) */
#define _ARF_POW_RSHIFT(r, a, len, s) \
    do { \
        slong __i; \
        ulong __s = (s); \
        FLINT_ASSERT(__s < 2 * FLINT_BITS); \
        if (__s >= FLINT_BITS) \
        { \
            for (__i = 0; __i < (len) - 1; __i++) \
                (r)[__i] = (a)[__i + 1]; \
            (r)[(len) - 1] = 0; \
            __s -= FLINT_BITS; \
        } \
        else \
        { \
            for (__i = 0; __i < (len); __i++) \
                (r)[__i] = (a)[__i]; \
        } \
        if (__s != 0) \
        { \
            for (__i = 0; __i < (len) - 1; __i++) \
                (r)[__i] = ((r)[__i] >> __s) | ((r)[__i + 1] << (FLINT_BITS - __s)); \
            (r)[(len) - 1] >>= __s; \
        } \
    } while (0)

FLINT_FORCE_INLINE int
_arf_small_pow_1(ulong * y2, ulong * y1, ulong * y0, slong * kout, ulong * err, ulong dx, slong ex, ulong dy, slong ey, int ysgn)
{
    ulong l3, l2, l1, l0, h, l, q, t2, t1, t0, eprop;
    ulong P[5], Z[5];
    slong k;
    int neg;

    if (ey > FLINT_BITS - 1 || ey < -(FLINT_BITS - 1))
        return ARF_ELEM_UNCERTAIN;

    if (ey <= 40)
    {
        /* Lambda at two fraction limbs suffices: (l3, l2, l1, 0) */
        neg = _arf_small_log_2(&l3, &l2, &l1, dx, 0, ex);
        l0 = 0;
    }
    else
        neg = _arf_small_log_3(&l3, &l2, &l1, &l0, dx, 0, 0, ex);

    /* P = Lambda my exactly, the value P / B^4 */
    umul_ppmm(P[1], P[0], l0, dy);
    umul_ppmm(h, l, l1, dy);
    add_ssaaaa(P[2], P[1], h, l, UWORD(0), P[1]);
    umul_ppmm(h, l, l2, dy);
    add_ssaaaa(P[3], P[2], h, l, UWORD(0), P[2]);
    umul_ppmm(h, l, l3, dy);
    add_ssaaaa(P[4], P[3], h, l, UWORD(0), P[3]);

    /* Z = floor(|t| B^3) = P >> (FLINT_BITS - ey) */
    _ARF_POW_RSHIFT(Z, P, 5, FLINT_BITS - ey);
    if (Z[4] != 0 || Z[3] >= (UWORD(1) << 61))
        return ARF_ELEM_UNCERTAIN;

    q = _arf_small_exp_reduce_3(&t2, &t1, &t0, Z[3], Z[2], Z[1], Z[0]);
    if (neg ^ ysgn)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 3);
        sub_dddmmmsss(t2, t1, t0, L[2], L[1], L[0], t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_2(y2, y1, y0, t2, t1);

    /* 2 (18 2^ey 2^64 + 1) 2^-64 units with Lambda at two limbs,
       2 (26 2^ey + 1) 2^-64 units at three */
    if (ey <= 40)
        eprop = (ey >= 0) ? (UWORD(36) << ey) + 2 : 38;
    else
        eprop = (UWORD(52) >> (FLINT_BITS - ey)) + 3;   /* (the shift floors) */
    *kout = k;
    *err = ARF_SMALL_EXP_ERR_2 + eprop;
    return 0;
}

FLINT_FORCE_INLINE int
_arf_small_pow_2(ulong * y3, ulong * y2, ulong * y1, ulong * y0, slong * kout, ulong * err,
    ulong dx1, ulong dx0, slong ex, ulong dy1, ulong dy0, slong ey, int ysgn)
{
    ulong l3, l2, l1, l0, h, l, q, t3, t2, t1, t0, eprop;
    ulong a4, a3, a2, a1, a0, b4, b3, b2, b1, b0;
    ulong P[6], Z[6];
    slong k;
    int neg;

    if (ey > 40 || ey < -(FLINT_BITS - 1))
        return ARF_ELEM_UNCERTAIN;

    neg = _arf_small_log_3(&l3, &l2, &l1, &l0, dx1, dx0, 0, ex);

    /* P = Lambda (my1, my0) exactly, the value P / B^5 */
    umul_ppmm(a1, a0, l0, dy0);
    umul_ppmm(h, l, l1, dy0);
    add_ssaaaa(a2, a1, h, l, UWORD(0), a1);
    umul_ppmm(h, l, l2, dy0);
    add_ssaaaa(a3, a2, h, l, UWORD(0), a2);
    umul_ppmm(h, l, l3, dy0);
    add_ssaaaa(a4, a3, h, l, UWORD(0), a3);

    umul_ppmm(b1, b0, l0, dy1);
    umul_ppmm(h, l, l1, dy1);
    add_ssaaaa(b2, b1, h, l, UWORD(0), b1);
    umul_ppmm(h, l, l2, dy1);
    add_ssaaaa(b3, b2, h, l, UWORD(0), b2);
    umul_ppmm(h, l, l3, dy1);
    add_ssaaaa(b4, b3, h, l, UWORD(0), b3);

    P[0] = a0;
    add_sssssaaaaaaaaaa(P[5], P[4], P[3], P[2], P[1], UWORD(0), a4, a3, a2, a1,
        b4, b3, b2, b1, b0);

    /* Z = floor(|t| B^4) = P >> (FLINT_BITS - ey) */
    _ARF_POW_RSHIFT(Z, P, 6, FLINT_BITS - ey);
    if (Z[5] != 0 || Z[4] >= (UWORD(1) << 61))
        return ARF_ELEM_UNCERTAIN;

    q = _arf_small_exp_reduce_4(&t3, &t2, &t1, &t0, Z[4], Z[3], Z[2], Z[1], Z[0]);
    if (neg ^ ysgn)
    {
        nn_srcptr L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, 4);
        sub_ddddmmmmssss(t3, t2, t1, t0, L[3], L[2], L[1], L[0], t3, t2, t1, t0);
        k = -(slong) (q + 1);
    }
    else
        k = (slong) q;

    _mp_real_small_exp_3(y3, y2, y1, y0, t3, t2, t1);

    /* 2 (26 2^ey + 1) units of B^-3 */
    eprop = (ey >= 0) ? (UWORD(52) << ey) + 2 : 54;
    *kout = k;
    *err = ARF_SMALL_EXP_ERR_3 + eprop;
    return 0;
}
#endif

#endif

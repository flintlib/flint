/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARF_MPN_SMALL_H
#define ARF_MPN_SMALL_H

#include <string.h>
#include "arf.h"
#include "mpn_extras/fp_small.h"

/*
    The quotient and root kernels of arf at one and two limbs, inline,
    for callers with fixed small sizes (nfloat, and the register paths of
    the kernels themselves): the register kernels of mpn_extras/fp_small.h
    with the rounding at a limb boundary (prec = rn FLINT_BITS) done on
    the exact floor and the remainder, or on a certified guard limb. The
    operands are single limbs resp. limb pairs (a1, a0), the fractions
    normalised (top bit set); z receives rn limbs and *zexp the correction
    to the natural exponent as for the kernels of arf.h. The _approx
    variants follow the contract of the _approx kernels.
*/

/* {z, rn} the exact floor (top bit set) of the value, whose bits below
   are summarised by half (the next bit) and sticky (whether any bit
   below that is set): rounded in place at prec = rn FLINT_BITS in the
   mode rnd, *zexp = 0 or 1 (the carry); returns the inexact flag */
FLINT_FORCE_INLINE int
_arf_mpn_finish_aligned(nn_ptr z, slong * zexp, slong rn, ulong half, ulong sticky, arf_rnd_t rnd)
{
    int up;

    if (rnd == ARF_RND_DOWN)
        up = 0;
    else if (rnd == ARF_RND_UP)
        up = (half | sticky) != 0;
    else
        up = half && (sticky || (z[0] & 1));

    *zexp = 0;
    if (up)
    {
        if (rn == 1)
        {
            z[0]++;
            if (z[0] == 0)
            {
                z[0] = LIMB_TOP;
                *zexp = 1;
            }
        }
        else if (rn == 2)
        {
            add_ssaaaa(z[1], z[0], z[1], z[0], UWORD(0), UWORD(1));
            if (z[1] == 0)
            {
                z[1] = LIMB_TOP;
                *zexp = 1;
            }
        }
        else if (mpn_add_1(z, z, rn, 1))
        {
            z[rn - 1] = LIMB_TOP;
            *zexp = 1;
        }
    }

    return (half | sticky) != 0;
}

/* the limb below an exact quotient from its remainder r < d: the half
   bit [2 r >= d] on top, a set bit below for 2 r != d, r != 0 (so that
   the limb is nonzero iff the quotient is inexact) */
FLINT_FORCE_INLINE ulong
_arf_mpn_remainder_code(ulong r, ulong d)
{
    ulong h = d - r;
    if (r > h)
        return LIMB_TOP | 1;
    if (r == h)
        return LIMB_TOP;
    return (r != 0);
}

/* the certification of an approximate guard limb g within bound units:
   away from 0 and B - 1, and from B / 2 for rounding to nearest */
FLINT_FORCE_INLINE int
_arf_mpn_guard_certified(ulong g, ulong bound, arf_rnd_t rnd)
{
    if (g <= bound || g >= UWORD_MAX - bound)
        return 0;
    if (rnd == ARF_RND_NEAR && g - (LIMB_TOP - bound) <= 2 * bound)
        return 0;
    return 1;
}

/* copies of results: a few limbs in registers, longer ones by memcpy
   (the plain loop of flint_mpn_copyi costs several instructions per
   limb) */
FLINT_FORCE_INLINE void
_arf_mpn_copy_limbs(nn_ptr d, nn_srcptr s, slong n)
{
    if (n <= 8)
        _flint_mpn_copy_small(d, s, n);
    else
        memcpy(d, s, n * sizeof(ulong));
}

/* Quotients *****************************************************************/

/* Division by a single normalised limb d by a chain of 2/1 divisions:
   (q, nq) = floor(x B^(nq - xn + 1) / (2^s d)) exactly for xn <= nq + 1
   (the top bit set; for xn = nq + 1 the bit shifted out of x is folded
   into the code), s = [x[xn-1] >= d]: x B^(nq - xn + 1) is divided
   first (nq + 1 quotient limbs, the top one 0 or 1) and the quotient
   shifted down by one bit when the top limb is set. Returns the limb
   below the quotient (its top bit the half bit, the rest the sticky
   bit, as _arf_mpn_remainder_code); q has nq + 1 limbs, t is scratch
   of nq + 1 limbs. Up to eight limbs flint_mpn_divrem_1_preinv is used
   with the inverse dinv = FLINT_MPN_INVERT_LIMB(d), computed here when
   0 is passed; beyond, GMP's mpn_divrem_1. */
FLINT_FORCE_INLINE ulong
_arf_mpn_div_chain(nn_ptr q, nn_ptr t, ulong * s, nn_srcptr x, slong xn, ulong d, ulong dinv, slong nq)
{
    slong nn = nq + 1;
    ulong r;

    if (nn <= 8)
    {
        if (dinv == 0)
            FLINT_MPN_INVERT_LIMB(dinv, d);
        _flint_mpn_zero_small(t, nn - xn);
        _flint_mpn_copy_small(t + nn - xn, x, xn);
        r = flint_mpn_divrem_1_preinv(q, t, nn, d, dinv, 0);
    }
    else
        r = mpn_divrem_1(q, nn - xn, x, xn, d);

    if (q[nq] != 0)
    {
        /* the quotient in [1, 2): halved, the bit shifted out as the
           half bit, the remainder below it */
        ulong bit = q[0] & 1;
        mpn_rshift(q, q, nn, 1);
        *s = 1;
        return (bit << (FLINT_BITS - 1)) | (r != 0);
    }

    *s = 0;
    return _arf_mpn_remainder_code(r, d);
}

/* x / d at n limbs (prec = n FLINT_BITS) by the chain, rounded resp.
   truncated; scratch of 2 n + 2 limbs */
FLINT_FORCE_INLINE int
_arf_mpn_div_limb(nn_ptr z, slong * zexp, nn_srcptr x, slong n, ulong d, ulong dinv, arf_rnd_t rnd, nn_ptr scratch)
{
    ulong code, s;
    int inexact;

    code = _arf_mpn_div_chain(scratch, scratch + n + 1, &s, x, n, d, dinv, n);
    _arf_mpn_copy_limbs(z, scratch, n);
    inexact = _arf_mpn_finish_aligned(z, zexp, n, code >> (FLINT_BITS - 1), code & 1, rnd);
    *zexp += s;
    return inexact;
}

FLINT_FORCE_INLINE void
_arf_mpn_div_limb_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong n, ulong d, ulong dinv, nn_ptr scratch)
{
    ulong s;

    _arf_mpn_div_chain(scratch, scratch + n + 1, &s, x, n, d, dinv, n);
    _arf_mpn_copy_limbs(z, scratch, n);
    *zexp = s;
}

/* x / y at one limb: one 2/1 division, exact with the remainder */
FLINT_FORCE_INLINE int
_arf_mpn_div_1(nn_ptr z, slong * zexp, ulong x, ulong y, arf_rnd_t rnd)
{
    ulong r, s, code;
    int inexact;
    z[0] = _flint_mpn_fdiv_1_kernel(&r, &s, x, y);
    code = _arf_mpn_remainder_code(r, y);
    inexact = _arf_mpn_finish_aligned(z, zexp, 1, code >> (FLINT_BITS - 1), code & 1, rnd);
    *zexp += s;
    return inexact;
}

FLINT_FORCE_INLINE void
_arf_mpn_div_1_approx(nn_ptr z, slong * zexp, ulong x, ulong y)
{
    ulong r, s;
    z[0] = _flint_mpn_fdiv_1_kernel(&r, &s, x, y);
    *zexp = s;
}

/* x / y at two limbs: two 3/2 divisions, exact with the inexact flag,
   which decides the directed modes; rounding to nearest needs the half
   bit, from a third quotient limb */
FLINT_FORCE_INLINE int
_arf_mpn_div_2(nn_ptr z, slong * zexp, ulong x1, ulong x0, ulong y1, ulong y0, arf_rnd_t rnd)
{
    ulong s, g;
    int inexact;

    if (rnd != ARF_RND_NEAR)
    {
        inexact = _flint_mpn_fdiv_2_kernel(z + 1, z, &s, x1, x0, y1, y0);
        _arf_mpn_finish_aligned(z, zexp, 2, 0, inexact, rnd);
    }
    else
    {
        inexact = _flint_mpn_fdiv_3_kernel(z + 1, z, &g, &s, x1, x0, y1, y0, 0);
        inexact = _arf_mpn_finish_aligned(z, zexp, 2, g >> (FLINT_BITS - 1), ((g << 1) != 0) | inexact, rnd);
    }
    *zexp += s;
    return inexact;
}

FLINT_FORCE_INLINE void
_arf_mpn_div_2_approx(nn_ptr z, slong * zexp, ulong x1, ulong x0, ulong y1, ulong y0)
{
    ulong s;
    _flint_mpn_fdiv_2_kernel(z + 1, z, &s, x1, x0, y1, y0);
    *zexp = s;
}

#if FLINT_BITS == 64

/* Square roots **************************************************************/

/* sqrt(x 2^-t) at one limb: the exact floor of the root of the two-limb
   radicand x B / 2^t with the remainder (the root exceeds the floor s by
   at least 1/2 iff rem > s: no ties) */
FLINT_FORCE_INLINE int
_arf_mpn_sqrt_1(nn_ptr z, slong * zexp, ulong x, int t, arf_rnd_t rnd)
{
    ulong r1, r0, half;
    int inexact;

    inexact = _flint_mpn_isqrt_2_rem(z, &r1, &r0, x >> t, (x << (FLINT_BITS - 1)) & (-(ulong) t));
    half = (r1 != 0 || r0 > z[0]);
    return _arf_mpn_finish_aligned(z, zexp, 1, half, inexact, rnd);
}

FLINT_FORCE_INLINE void
_arf_mpn_sqrt_1_approx(nn_ptr z, slong * zexp, ulong x, int t)
{
    _flint_mpn_isqrt_2(z, x >> t, (x << (FLINT_BITS - 1)) & (-(ulong) t));
    *zexp = 0;
}

/* sqrt(x 2^-t) at two limbs: three limbs within FLINT_MPN_SMALL_SQRT_ERR_3
   units, certified when the guard limb lies away from the rounding
   boundaries, else the exact root */
FLINT_FORCE_INLINE int
_arf_mpn_sqrt_2(nn_ptr z, slong * zexp, ulong x1, ulong x0, int t, arf_rnd_t rnd)
{
    ulong s2, s1, s0, x[2];
    int code;

    _flint_mpn_small_sqrt_3(&s2, &s1, &s0, x1, x0, 0, t);
    if (_arf_mpn_guard_certified(s0, FLINT_MPN_SMALL_SQRT_ERR_3, rnd))
    {
        z[0] = s1;
        z[1] = s2;
        _arf_mpn_finish_aligned(z, zexp, 2, s0 >> (FLINT_BITS - 1), 1, rnd);
        return 1;
    }

    x[0] = x0;
    x[1] = x1;
    if (rnd == ARF_RND_NEAR)
        code = flint_mpn_sqrthigh_exact_half(z, x, 2, t, 1);
    else
        code = flint_mpn_sqrthigh_exact(z, x, 2, t, 1);
    return _arf_mpn_finish_aligned(z, zexp, 2, code >> 1, code & 1, rnd);
}

FLINT_FORCE_INLINE void
_arf_mpn_sqrt_2_approx(nn_ptr z, slong * zexp, ulong x1, ulong x0, int t)
{
    ulong s2, s1, s0, x[2];

    *zexp = 0;
    _flint_mpn_small_sqrt_3(&s2, &s1, &s0, x1, x0, 0, t);
    if (s0 > FLINT_MPN_SMALL_SQRT_ERR_3 && s0 < UWORD_MAX - FLINT_MPN_SMALL_SQRT_ERR_3)
    {
        z[0] = s1;
        z[1] = s2;
        return;
    }

    x[0] = x0;
    x[1] = x1;
    flint_mpn_sqrthigh_exact(z, x, 2, t, 1);
}

/* Reciprocal square roots ***************************************************/

/* the exact fallback: the floor at rn + 1 limbs (irrational: inexact)
   rounded; x must not be a power of two with t = 1 (M = 1), which the
   callers handle (1/2 at the exponent one up) */
FLINT_FORCE_INLINE int
_arf_mpn_rsqrt_small_exact(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong rn, arf_rnd_t rnd)
{
    ulong q[3];
    flint_mpn_rsqrthigh_exact(q, x, xn, t, rn);
    if (rn == 1)
        z[0] = q[1];
    else
    {
        z[0] = q[1];
        z[1] = q[2];
    }
    _arf_mpn_finish_aligned(z, zexp, rn, q[0] >> (FLINT_BITS - 1), 1, rnd);
    return 1;
}

/* 1 / sqrt(x 2^-t) at one limb: two limbs within 6 units (the third-order
   refinement), certified away from the rounding boundaries; a power of
   two with t = 1 is the caller's (1/2 at the exponent one up) */
FLINT_FORCE_INLINE int
_arf_mpn_rsqrt_1(nn_ptr z, slong * zexp, ulong x, int t, arf_rnd_t rnd)
{
    ulong y1, y0;

    _flint_mpn_rsqrt_1_approx_128(&y1, &y0, x, t, 1);
    if (_arf_mpn_guard_certified(y0, 8, rnd))
    {
        z[0] = y1;
        _arf_mpn_finish_aligned(z, zexp, 1, y0 >> (FLINT_BITS - 1), 1, rnd);
        return 1;
    }
    return _arf_mpn_rsqrt_small_exact(z, zexp, &x, 1, t, 1, rnd);
}

FLINT_FORCE_INLINE void
_arf_mpn_rsqrt_1_approx(nn_ptr z, slong * zexp, ulong x, int t)
{
    _flint_mpn_rsqrt_1_approx(z, x, t, 0);
    *zexp = 0;
}

/* 1 / sqrt(x 2^-t) at two limbs: three limbs within
   FLINT_MPN_SMALL_RSQRT_ERR_3 units */
FLINT_FORCE_INLINE int
_arf_mpn_rsqrt_2(nn_ptr z, slong * zexp, ulong x1, ulong x0, int t, arf_rnd_t rnd)
{
    ulong z2, z1, z0, x[2];

    _flint_mpn_small_rsqrt_3(&z2, &z1, &z0, x1, x0, t);
    if (_arf_mpn_guard_certified(z0, FLINT_MPN_SMALL_RSQRT_ERR_3, rnd))
    {
        z[0] = z1;
        z[1] = z2;
        _arf_mpn_finish_aligned(z, zexp, 2, z0 >> (FLINT_BITS - 1), 1, rnd);
        return 1;
    }
    x[0] = x0;
    x[1] = x1;
    return _arf_mpn_rsqrt_small_exact(z, zexp, x, 2, t, 2, rnd);
}

FLINT_FORCE_INLINE void
_arf_mpn_rsqrt_2_approx(nn_ptr z, slong * zexp, ulong x1, ulong x0, int t)
{
    ulong z2, z1, z0;
    _flint_mpn_small_rsqrt_3(&z2, &z1, &z0, x1, x0, t);
    z[0] = z1;
    z[1] = z2;
    *zexp = 0;
}

#endif

#endif

/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"
#include "arf/mpn_small.h"

/*
    Quotients, reciprocals, square roots and reciprocal square roots on
    limb-array mantissas (the _arf_mpn_ interface of arf.h), in two
    flavours: rounded in a sign-free mode (ARF_RND_DOWN, ARF_RND_UP or
    ARF_RND_NEAR) with the inexact flag returned, and _approx (rn limbs
    within an ulp of the magnitude, nothing else: the kernels of nfloat
    without directed rounding).

    Small sizes go through the register kernels of mpn_extras/fp_small.h
    (one and two result limbs, inline in arf/mpn_small.h for callers
    with fixed sizes; a divisor of a single limb, such as an integer,
    through a chain of 2/1 divisions at any size), which give exact
    floors with the remainder (quotients, the one-limb root) or a
    certified approximation (the roots at two limbs). Longer quotients
    use the guard-limb kernels of mpn_extras (flint_mpn_divhigh and
    friends): n + 1 limbs with the top bit set, the lowest a guard limb
    within a small bound of the exact value in either direction, and the
    exact floor with an inexact flag as the fallback; the _approx
    quotient uses the schoolbook approximate division at rn limbs
    directly. Square roots are exact integer roots of the mantissa padded
    to 2 rn limbs (flint_mpn_sqrthigh_exact, with the half bit for
    rounding to nearest at a limb boundary), which is cheaper than the
    guard limb up to the Newton range; long reciprocal roots keep the
    guard limb.

    The rounded kernels certify the rounding from the guard window
    (_arf_mpn_round_guarded) and fall back to the exact floor with a
    sticky bit, so that an exact quotient or root is reported exact. The
    kernels write rn = ceil(prec / FLINT_BITS) limbs into z, which may
    alias an operand, and use a buffer of rn + 1 limbs of their own.
*/

/* the exponents of the roots: x = A 2^e with A in [1/2, 1), t = e mod 2;
   sqrt(x) = S 2^((e + t) / 2) and 1/sqrt(x) = M 2^((1 - t) - (e - t) / 2)
   for the S and M of flint_mpn_sqrthigh / rsqrthigh: the callers hold
   those exponents, the kernels return the corrections */

/* zero low limbs of a mantissa (nfloat pads to a fixed length) stripped:
   the register and chain paths see the true length */
#define STRIP_LOW(d, n) do { while ((d)[0] == 0) { (d)++; (n)--; } } while (0)

/* whether the normalised mantissa is a power of two (1/2): the top limb
   alone set */
FLINT_FORCE_INLINE int
_is_pow2(nn_srcptr d, slong n)
{
    return d[n - 1] == LIMB_TOP && (n == 1 || flint_mpn_zero_p(d, n - 1));
}

/* {x, xn} (normalised) rounded to prec bits into {z, rn}, for a quotient
   by a power of two or a root that is one (exact when it fits: zero
   padded) */
static int
_round_limbs(nn_ptr z, slong * zexp, nn_srcptr xd, slong xn, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (xn * FLINT_BITS <= prec)
    {
        /* (z may alias x: the copy in the safe direction) */
        if (z + rn - xn > xd)
            flint_mpn_copyd(z + rn - xn, xd, xn);
        else
            flint_mpn_copyi(z + rn - xn, xd, xn);
        _flint_mpn_zero_small(z, rn - xn);
        *zexp = 0;
        return 0;
    }

    return _arf_mpn_set_round(z, zexp, xd, xn, prec, rnd);
}

/* the power of two 1/2 into {z, rn} */
static void
_set_half(nn_ptr z, slong rn)
{
    _flint_mpn_zero_small(z, rn - 1);
    z[rn - 1] = LIMB_TOP;
}

/* The exact floor {q, qn} (top bit set) of a value with prec <= qn
   FLINT_BITS, the limb q[-1] encoding the discarded part below q: its
   top bit the half bit and its low bit set iff there are further bits
   (so that q[-1] != 0 iff the floor is inexact; a sticky bit alone,
   without the half bit, suffices when the rounding position lies inside
   q). Rounds to prec bits into {z, rn}. */
FLINT_FORCE_INLINE int
_round_floor(nn_ptr z, slong * zexp, nn_ptr q, slong qn, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    if (prec == rn * FLINT_BITS)
    {
        slong gn = qn - rn, i;
        ulong g = q[gn - 1], sticky;

        sticky = (g << 1) != 0;
        for (i = gn - 2; i >= -1; i--)
            sticky |= (q[i] != 0);
        _arf_mpn_copy_limbs(z, q + gn, rn);
        return _arf_mpn_finish_aligned(z, zexp, rn, g >> (FLINT_BITS - 1), sticky, rnd);
    }

    return _arf_mpn_set_round(z, zexp, q - 1, qn + 1, prec, rnd);
}

/* Quotients *****************************************************************/

/* [A >= Y] for fractions of an resp. bn limbs (zero-extended) */
FLINT_FORCE_INLINE ulong
_fraction_ge(nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    slong m;
    int c;

    if (a[an - 1] != b[bn - 1])
        return a[an - 1] > b[bn - 1];

    m = FLINT_MIN(an, bn);
    c = mpn_cmp(a + an - m, b + bn - m, m);
    if (c != 0)
        return c > 0;
    if (bn > an)
        return flint_mpn_zero_p(b, bn - an);
    return 1;
}

/* (t, nn) = floor(x B^(nn - xn) / 2^s) for s in {0, 1} (the bit shifted
   out lands in the padding, lost when nn = xn), or the top nn limbs of
   x / 2^s when xn > nn */
FLINT_FORCE_INLINE void
_shifted_numerator(nn_ptr t, slong nn, nn_srcptr x, slong xn, ulong s)
{
    slong off, i;
    ulong mask = -s, lo = 0;

    if (xn > nn)
    {
        x += xn - nn;
        xn = nn;
    }
    off = nn - xn;
    _flint_mpn_zero_small(t, off);

    for (i = xn - 1; i >= 0; i--)
    {
        t[off + i] = (x[i] >> s) | lo;
        lo = (x[i] << (FLINT_BITS - 1)) & mask;
    }
    if (off > 0)
        t[off - 1] = lo;
}

/* the register paths of the quotients: result limbs rn <= 2 with xn <= 2
   and yn <= 3 (exact: the floor with the remainder) */
#define DIV_REGISTER(rn, xn, yn) ((rn) <= 2 && (xn) <= 2 && (yn) <= 3)

/* The exact floor (q, nq) = floor(x B^(nq + yn - xn) / (2^s y)) with its
   remainder by the register-based division (_flint_mpn_tdiv_qr_small)
   for divisors of up to three limbs (where flint_mpn_divhigh would use
   it too, at one more limb and with the certification to follow),
   s = [X >= Y]; requires
   xn < nq + yn (the bit shifted out of x lands in the padding). Returns
   the limb below the quotient (the half bit on top, a low bit for
   further bits), as _arf_mpn_remainder_code: the remainder R against
   Y - R. q has nq + 1 limbs (the top one comes out zero), t is scratch
   of nq + 3 yn limbs (the numerator, the remainder and Y - R). */
#define DIV_SMALL_EXACT(nq, xn, yn) ((yn) <= 3 && (xn) < (nq) + (yn))

FLINT_FORCE_INLINE ulong
_div_small_exact(nn_ptr q, nn_ptr t, ulong * s, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong nq)
{
    slong nn = nq + yn;
    nn_ptr r = t + nn, d = r + yn;
    int c;

    *s = _fraction_ge(x, xn, y, yn);
    _shifted_numerator(t, nn, x, xn, *s);
    _flint_mpn_tdiv_qr_small(q, r, t, nn, y, yn);
    FLINT_ASSERT(q[nq] == 0 && (q[nq - 1] >> (FLINT_BITS - 1)));

    if (flint_mpn_zero_p(r, yn))
        return 0;
    /* 2 R against Y: R against Y - R */
    mpn_sub_n(d, y, r, yn);
    c = mpn_cmp(r, d, yn);
    if (c > 0)
        return LIMB_TOP | 1;
    if (c == 0)
        return LIMB_TOP;
    return 1;
}

/* the schoolbook approximate division of the _approx quotient up to
   this many result limbs (as DIVHIGH_BASECASE_MAX of divhigh.c) */
#define DIV_APPROX_BASECASE_MAX 80

/* (q, 2 or 3) = floor(x B^qn / (2^s y)) for xn <= 2 and yn <= 3 in
   registers, qn = 3 limbs when yn = 3 or three limbs are requested;
   q[-1] receives the inexact flag */
FLINT_FORCE_INLINE slong
_div_register(nn_ptr q, ulong * s, nn_srcptr x, slong xn, nn_srcptr y, slong yn, int three)
{
    ulong x1 = x[xn - 1], x0 = (xn == 2) ? x[0] : 0;
    ulong y1 = y[yn - 1], y0 = (yn >= 2) ? y[yn - 2] : 0;

    if (yn == 3 || three)
    {
        q[-1] = _flint_mpn_fdiv_3_kernel(q + 2, q + 1, q, s, x1, x0, y1, y0, (yn == 3) ? y[0] : 0);
        return 3;
    }
    else
    {
        q[-1] = _flint_mpn_fdiv_2_kernel(q + 1, q, s, x1, x0, y1, y0);
        return 2;
    }
}

int
_arf_mpn_div(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, adj;
    ulong s;
    nn_ptr q;
    int r;
    ARF_MPN_TMP_DECL

    FLINT_ASSERT(prec >= 1);

    STRIP_LOW(x, xn);
    STRIP_LOW(y, yn);

    if (yn == 1)
    {
        ulong d = y[0];

        if (d == LIMB_TOP)
        {
            /* x / 2^k */
            r = _round_limbs(z, zexp, x, xn, prec, rnd);
            *zexp += 1;
            return r;
        }

        if (rn == 1 && xn == 1)
        {
            if (prec == FLINT_BITS)
                return _arf_mpn_div_1(z, zexp, x[0], d, rnd);
            else
            {
                ulong qq[2], rem;
                qq[1] = _flint_mpn_fdiv_1_kernel(&rem, &s, x[0], d);
                qq[0] = _arf_mpn_remainder_code(rem, d);
                r = _round_floor(z, &adj, qq + 1, 1, prec, rnd);
                *zexp = (slong) s + adj;
                return r;
            }
        }

        if (xn <= rn)
        {
            /* the exact quotient at rn limbs with its remainder */
            ARF_MPN_TMP_ALLOC(q, 2 * rn + 3)
            q[0] = _arf_mpn_div_chain(q + 1, q + rn + 2, &s, x, xn, d, 0, rn);
            r = _round_floor(z, &adj, q + 1, rn, prec, rnd);
            *zexp = (slong) s + adj;
            return r;
        }
    }

    if (DIV_REGISTER(rn, xn, yn))
    {
        if (rn == 2 && xn == 2 && yn == 2 && prec == 2 * FLINT_BITS)
            return _arf_mpn_div_2(z, zexp, x[1], x[0], y[1], y[0], rnd);
        else
        {
            /* the exact floor at two limbs, or three when the half bit
               of a two-limb result is needed */
            ulong qq[4];
            slong qn = _div_register(qq + 1, &s, x, xn, y, yn,
                rn == 2 && rnd == ARF_RND_NEAR && prec == 2 * FLINT_BITS);
            r = _round_floor(z, &adj, qq + 1, qn, prec, rnd);
            *zexp = (slong) s + adj;
            return r;
        }
    }

    if (DIV_SMALL_EXACT(rn, xn, yn))
    {
        /* the exact floor at rn limbs with its remainder */
        ARF_MPN_TMP_ALLOC(q, 2 * rn + 3 * yn + 2)
        q[0] = _div_small_exact(q + 1, q + rn + 2, &s, x, xn, y, yn, rn);
        r = _round_floor(z, &adj, q + 1, rn, prec, rnd);
        *zexp = (slong) s + adj;
        return r;
    }

    ARF_MPN_TMP_ALLOC(q, rn + 1)
    s = flint_mpn_divhigh(q, x, xn, y, yn, rn);

    if (_arf_mpn_round_guarded(z, &adj, q + 1, rn, q[0], FLINT_BIT_COUNT(FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND), prec, rnd, 0) >= 0)
    {
        *zexp = (slong) s + adj;
        return 1;
    }

    r = flint_mpn_divhigh_exact(q, &s, x, xn, y, yn, rn);
    q[0] |= r;
    r = _arf_mpn_set_round(z, &adj, q, rn + 1, prec, rnd);
    *zexp = (slong) s + adj;
    return r;
}


void
_arf_mpn_div_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    ulong s;
    nn_ptr q;
    ARF_MPN_TMP_DECL

    FLINT_ASSERT(prec >= 1);

    STRIP_LOW(x, xn);
    STRIP_LOW(y, yn);

    if (yn == 1)
    {
        ulong d = y[0];

        if (d == LIMB_TOP)
        {
            if (xn >= rn)
                flint_mpn_copyi(z, x + xn - rn, rn);
            else
            {
                if (z + rn - xn > x)
                    flint_mpn_copyd(z + rn - xn, x, xn);
                else
                    flint_mpn_copyi(z + rn - xn, x, xn);
                _flint_mpn_zero_small(z, rn - xn);
            }
            *zexp = 1;
            return;
        }

        if (rn == 1 && xn == 1)
        {
            _arf_mpn_div_1_approx(z, zexp, x[0], d);
            return;
        }

        /* the chain on the top rn + 1 limbs of x (the rest below
           2 / B ulp of the quotient) */
        if (xn > rn + 1)
        {
            x += xn - (rn + 1);
            xn = rn + 1;
        }

        {
            ARF_MPN_TMP_ALLOC(q, 2 * rn + 2)
            _arf_mpn_div_chain(q, q + rn + 1, &s, x, xn, d, 0, rn);
            _arf_mpn_copy_limbs(z, q, rn);
            *zexp = (slong) s;
            return;
        }
    }

    if (DIV_REGISTER(rn, xn, yn))
    {
        if (rn == 2 && xn == 2 && yn == 2)
            _arf_mpn_div_2_approx(z, zexp, x[1], x[0], y[1], y[0]);
        else
        {
            ulong qq[4];
            slong qn = _div_register(qq + 1, &s, x, xn, y, yn, 0);
            /* (rn <= 2 here: the limbs copied directly, which also tells
               the compiler that the reads stay within qq) */
            z[0] = qq[1 + qn - rn];
            if (rn == 2)
                z[1] = qq[2 + qn - rn];
            *zexp = (slong) s;
        }
        return;
    }

    if (rn <= DIV_APPROX_BASECASE_MAX && yn >= 3)
    {
        /* the schoolbook approximate division at rn limbs: the floor
           or one more (the divisor truncated to rn + 2 limbs, the
           numerator x B^yn / 2^s to rn + yn limbs: both below an ulp) */
        slong nn;
        ulong dinv, top;
        nn_ptr t, qq;

        if (yn > rn + 2)
        {
            y += yn - (rn + 2);
            yn = rn + 2;
        }
        s = _fraction_ge(x, xn, y, yn);
        nn = rn + yn;
        ARF_MPN_TMP_ALLOC(t, nn + rn + 1)
        _shifted_numerator(t, nn, x, xn, s);
        dinv = flint_mpn_preinv1(y[yn - 1], y[yn - 2]);

        /* in place unless the output overlaps the divisor */
        qq = (z + rn > y && y + yn > z) ? t + nn : z;
        top = _flint_mpn_divapprox_basecase_preinv1(qq, t, nn, y, yn, dinv);

        /* the quotient lies in [1/2, 1) up to the truncations: clamped */
        if (FLINT_UNLIKELY(top != 0))
            flint_mpn_store(qq, rn, UWORD_MAX);
        else if (FLINT_UNLIKELY(!(qq[rn - 1] >> (FLINT_BITS - 1))))
            _set_half(qq, rn);
        if (qq != z)
            _arf_mpn_copy_limbs(z, qq, rn);
        *zexp = (slong) s;
        return;
    }

    /* the guard-limb quotient truncated (two-limb divisors: the exact
       register division) */
    ARF_MPN_TMP_ALLOC(q, rn + 1)
    s = flint_mpn_divhigh(q, x, xn, y, yn, rn);
    _arf_mpn_copy_limbs(z, q + 1, rn);
    *zexp = (slong) s;
}

/* Reciprocals ***************************************************************/

/* 1 = (1/2) 2^1: the quotient kernels with x = 1/2, the exponent one up */

int
_arf_mpn_inv(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec, arf_rnd_t rnd)
{
    ulong half = LIMB_TOP;
    int r = _arf_mpn_div(z, zexp, &half, 1, y, yn, prec, rnd);
    *zexp += 1;
    return r;
}


void
_arf_mpn_inv_approx(nn_ptr z, slong * zexp, nn_srcptr y, slong yn, slong prec)
{
    ulong half = LIMB_TOP;
    _arf_mpn_div_approx(z, zexp, &half, 1, y, yn, prec);
    *zexp += 1;
}

/* Square roots **************************************************************/

/* the exact root (sqrtrem of 2 rn limbs) is used where sqrthigh would use
   sqrtrem anyway, i.e. below the Newton range */
#define SQRT_EXACT_MAX (FLINT_MPN_SQRTREM_NEWTON_CUTOFF / 4)

/* ... but with the remainder (for the half bit) only up to this many
   limbs, beyond which the guard limb is cheaper (measured: 24-32 limbs
   on Skylake) */
#define SQRT_HALF_MAX 24

int
_arf_mpn_sqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, adj;
    nn_ptr q;
    int r;
    ARF_MPN_TMP_DECL

    FLINT_ASSERT(prec >= 1 && (t == 0 || t == 1));

    STRIP_LOW(x, xn);

    if (t && _is_pow2(x, xn))
    {
        /* x = 2^(e - 1) with e - 1 even: the root is a power of two */
        _set_half(z, rn);
        *zexp = 0;
        return 0;
    }

#if FLINT_BITS == 64
    if (rn == 1 && xn <= 2)
    {
        if (xn == 1 && prec == FLINT_BITS)
            return _arf_mpn_sqrt_1(z, zexp, x[0], t, rnd);
        else
        {
            /* the exact floor of the root of the two-limb radicand
               floor(A B^2 / 2^t) with the remainder: the root exceeds
               the floor s by at least 1/2 iff rem > s, or rem = s and
               the bit shifted out of A is set (no ties) */
            ulong x1 = x[xn - 1], x0 = (xn == 2) ? x[0] : 0, n1, n0, bit;
            ulong qq[2], r1, r0;
            int inexact;

            n1 = x1 >> t;
            n0 = (x0 >> t) | ((x1 << (FLINT_BITS - 1)) & (-(ulong) t));
            bit = x0 & t;
            inexact = _flint_mpn_isqrt_2_rem(qq + 1, &r1, &r0, n1, n0) | bit;
            if (r1 != 0 || r0 > qq[1] || (r0 == qq[1] && bit))
                qq[0] = LIMB_TOP | 1;
            else
                qq[0] = inexact;
            return _round_floor(z, zexp, qq + 1, 1, prec, rnd);
        }
    }

    if (rn == 2 && xn <= 2 && prec == 2 * FLINT_BITS)
        return _arf_mpn_sqrt_2(z, zexp, x[xn - 1], (xn == 2) ? x[0] : 0, t, rnd);

    if (rn == 2 && xn <= 3 && prec < 2 * FLINT_BITS)
    {
        /* three limbs within FLINT_MPN_SMALL_SQRT_ERR_3 units */
        ulong qq[2], s0;

        _flint_mpn_small_sqrt_3(qq + 1, qq, &s0, x[xn - 1], (xn >= 2) ? x[xn - 2] : 0, (xn == 3) ? x[0] : 0, t);
        if (_arf_mpn_round_guarded(z, zexp, qq, 2, s0, 2, prec, rnd, 0) >= 0)
            return 1;
        /* uncertified: the exact root below */
    }
#endif

    if (rn < SQRT_EXACT_MAX && (rnd != ARF_RND_NEAR || prec < rn * FLINT_BITS || rn <= SQRT_HALF_MAX))
    {
        /* the exact floor: with the half bit when rounding to nearest
           at a limb boundary, otherwise with the inexact flag alone */
        int code;

        if (prec == rn * FLINT_BITS)
        {
            if (rnd == ARF_RND_NEAR)
                code = flint_mpn_sqrthigh_exact_half(z, x, xn, t, rn - 1);
            else
                code = flint_mpn_sqrthigh_exact(z, x, xn, t, rn - 1);
            return _arf_mpn_finish_aligned(z, zexp, rn, code >> 1, code & 1, rnd);
        }

        ARF_MPN_TMP_ALLOC(q, rn + 1)
        code = flint_mpn_sqrthigh_exact(q + 1, x, xn, t, rn - 1);
        q[0] = code;
        return _arf_mpn_set_round(z, zexp, q, rn + 1, prec, rnd);
    }

    ARF_MPN_TMP_ALLOC(q, rn + 1)
    flint_mpn_sqrthigh(q, x, xn, t, rn);

    if (_arf_mpn_round_guarded(z, &adj, q + 1, rn, q[0], FLINT_BIT_COUNT(FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND), prec, rnd, 0) >= 0)
    {
        *zexp = adj;
        return 1;
    }

    r = flint_mpn_sqrthigh_exact(q, x, xn, t, rn);
    q[0] |= r;
    r = _arf_mpn_set_round(z, &adj, q, rn + 1, prec, rnd);
    *zexp = adj;
    return r;
}


void
_arf_mpn_sqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    FLINT_ASSERT(prec >= 1 && (t == 0 || t == 1));

    STRIP_LOW(x, xn);
    *zexp = 0;

#if FLINT_BITS == 64
    if (rn == 1)
    {
        /* the exact floor from the top two limbs (floor(sqrt(floor(N)))
           = floor(sqrt(N))) */
        ulong x1 = x[xn - 1], x0 = (xn >= 2) ? x[xn - 2] : 0;
        _flint_mpn_isqrt_2(z, x1 >> t, (x0 >> t) | ((x1 << (FLINT_BITS - 1)) & (-(ulong) t)));
        return;
    }

    if (rn == 2 && xn <= 2)
    {
        _arf_mpn_sqrt_2_approx(z, zexp, x[xn - 1], (xn == 2) ? x[0] : 0, t);
        return;
    }
#endif

    if (rn < SQRT_EXACT_MAX)
    {
        flint_mpn_sqrthigh_exact(z, x, xn, t, rn - 1);
        return;
    }

    /* the guard-limb root truncated */
    {
        nn_ptr q;
        ARF_MPN_TMP_DECL
        ARF_MPN_TMP_ALLOC(q, rn + 1)
        flint_mpn_sqrthigh(q, x, xn, t, rn);
        _arf_mpn_copy_limbs(z, q + 1, rn);
    }
}

/* Reciprocal square roots ***************************************************/

/* the guard-limb root: the Newton iteration directly in its range (the
   dispatch of flint_mpn_rsqrthigh skipped) */
FLINT_FORCE_INLINE void
_rsqrthigh(nn_ptr q, nn_srcptr x, slong xn, int t, slong rn)
{
    if (rn >= 2 && rn + 1 <= FLINT_MPN_RSQRTHIGH_MP_REAL_CUTOFF)
    {
        _flint_mpn_rsqrthigh_newton(q, rn + 1, x, xn, t);
        /* the top bit is set unless the approximation dipped below 1/2 */
        if (FLINT_UNLIKELY(!(q[rn] >> (FLINT_BITS - 1))))
            _set_half(q, rn + 1);
    }
    else
        flint_mpn_rsqrthigh(q, x, xn, t, rn);
}

int
_arf_mpn_rsqrt(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec, arf_rnd_t rnd)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS, adj;
    nn_ptr q;
    ARF_MPN_TMP_DECL

    FLINT_ASSERT(prec >= 1 && (t == 0 || t == 1));

    STRIP_LOW(x, xn);

    if (t && _is_pow2(x, xn))
    {
        /* x = 2^(e - 1) with e - 1 even: 1/sqrt(x) = 2^(-(e - 1)/2), i.e.
           1/2 at the exponent 1 - (e - 1)/2 = RSQRT_EXP + 1 */
        _set_half(z, rn);
        *zexp = 1;
        return 0;
    }

#if FLINT_BITS == 64
    if (rn == 1 && xn == 1)
    {
        if (prec == FLINT_BITS)
            return _arf_mpn_rsqrt_1(z, zexp, x[0], t, rnd);
        else
        {
            /* two limbs within 6 units (the third-order refinement) */
            ulong y1, y0;
            _flint_mpn_rsqrt_1_approx_128(&y1, &y0, x[0], t, 1);
            if (_arf_mpn_round_guarded(z, zexp, &y1, 1, y0, 3, prec, rnd, 0) >= 0)
                return 1;
        }
    }
    else if (rn == 2 && xn <= 2)
    {
        if (prec == 2 * FLINT_BITS)
            return _arf_mpn_rsqrt_2(z, zexp, x[xn - 1], (xn == 2) ? x[0] : 0, t, rnd);
        else
        {
            /* three limbs within FLINT_MPN_SMALL_RSQRT_ERR_3 units */
            ulong qq[2], z0;
            _flint_mpn_small_rsqrt_3(qq + 1, qq, &z0, x[xn - 1], (xn == 2) ? x[0] : 0, t);
            if (_arf_mpn_round_guarded(z, zexp, qq, 2, z0, 2, prec, rnd, 0) >= 0)
                return 1;
        }
    }
    else
#endif
    {
        ARF_MPN_TMP_ALLOC(q, rn + 1)
        _rsqrthigh(q, x, xn, t, rn);

        if (_arf_mpn_round_guarded(z, &adj, q + 1, rn, q[0], FLINT_BIT_COUNT(FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND), prec, rnd, 0) >= 0)
        {
            *zexp = adj;
            return 1;
        }
    }

    /* uncertified (the reciprocal root of a number that is not a power
       of two is irrational: inexact) */
    ARF_MPN_TMP_ALLOC(q, rn + 1)
    flint_mpn_rsqrthigh_exact(q, x, xn, t, rn);
    q[0] |= 1;
    _arf_mpn_set_round(z, &adj, q, rn + 1, prec, rnd);
    *zexp = adj;
    return 1;
}


void
_arf_mpn_rsqrt_approx(nn_ptr z, slong * zexp, nn_srcptr x, slong xn, int t, slong prec)
{
    slong rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    nn_ptr q;
    ARF_MPN_TMP_DECL

    FLINT_ASSERT(prec >= 1 && (t == 0 || t == 1));

    STRIP_LOW(x, xn);

    if (t && _is_pow2(x, xn))
    {
        _set_half(z, rn);
        *zexp = 1;
        return;
    }

#if FLINT_BITS == 64
    if (rn == 1 && xn == 1)
    {
        _arf_mpn_rsqrt_1_approx(z, zexp, x[0], t);
        return;
    }

    if (rn == 2 && xn <= 2)
    {
        _arf_mpn_rsqrt_2_approx(z, zexp, x[xn - 1], (xn == 2) ? x[0] : 0, t);
        return;
    }
#endif

    /* the guard-limb root truncated */
    *zexp = 0;
    ARF_MPN_TMP_ALLOC(q, rn + 1)
    _rsqrthigh(q, x, xn, t, rn);
    _arf_mpn_copy_limbs(z, q + 1, rn);
}

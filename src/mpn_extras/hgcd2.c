/*
    From GMP's mpn/generic/hgcd2.c:
    Copyright 1996, 1998, 2000-2004, 2008, 2012, 2019, 2020 Free Software
    Foundation, Inc.

    From GMP's mpn/generic/matrix22_mul1_inverse_vector.c:
    Copyright 2008, 2010 Free Software Foundation, Inc.
    Contributed to the GNU project by Niels Möller.

    From GMP's mpn/generic/gcd_22.c:
    Copyright 2019 Free Software Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Double-limb half gcd (hgcd2) and the 2 x 2 single-limb matrix-vector
    products of Lehmer's algorithm, ported from GMP (the products with a
    different, single-pass loop; see below). The division helpers are
    GMP's methods 3 and 5 for div1 (selected by
    FLINT_MPN_HGCD2_DIV1_METHOD) and method 2 for div2.
*/

/* Method 3 measured about 7% faster than method 5 on Skylake-class Xeons
   (contrary to GMP's tuning), so it is the default. */
#ifndef FLINT_MPN_HGCD2_DIV1_METHOD
# define FLINT_MPN_HGCD2_DIV1_METHOD 3
#endif

typedef struct
{
    mp_limb_t d0, d1;
}
_double_limb_t;

#if FLINT_MPN_HGCD2_DIV1_METHOD == 3

/* Single-limb division optimized for small quotients: d1 = q, d0 = r. */
FLINT_FORCE_INLINE _double_limb_t
div1(mp_limb_t n0, mp_limb_t d0)
{
    _double_limb_t res;

    if (FLINT_UNLIKELY((d0 >> (FLINT_BITS - 3)) != 0)
        || FLINT_UNLIKELY(n0 >= (d0 << 3)))
    {
        res.d1 = n0 / d0;
        res.d0 = n0 - res.d1 * d0;
    }
    else
    {
        mp_limb_t q, mask;

        d0 <<= 2;

        mask = -(mp_limb_t) (n0 >= d0);
        n0 -= d0 & mask;
        q = 4 & mask;

        d0 >>= 1;
        mask = -(mp_limb_t) (n0 >= d0);
        n0 -= d0 & mask;
        q += 2 & mask;

        d0 >>= 1;
        mask = -(mp_limb_t) (n0 >= d0);
        n0 -= d0 & mask;
        q -= mask;

        res.d0 = n0;
        res.d1 = q;
    }
    return res;
}

#elif FLINT_MPN_HGCD2_DIV1_METHOD == 5

/* Table of inverses of the 7 most significant divisor bits; full division
   is needed about 0.49% of the time. */
static const unsigned char div1_tab[64] = {
    255,251,247,243,239,236,233,229,226,223,220,217,214,211,209,206,
    203,201,198,196,194,191,189,187,185,183,181,179,177,175,173,171,
    169,167,166,164,162,161,159,158,156,155,153,152,150,149,147,146,
    145,143,142,141,140,139,137,136,135,134,133,132,131,130,129,128
};

#define DIV1_NBITS 7

FLINT_FORCE_INLINE _double_limb_t
div1(mp_limb_t n0, mp_limb_t d0)
{
    int ncnt, dcnt;
    size_t dbi;
    mp_limb_t inv, q0, r0, mask;
    _double_limb_t res;

    ncnt = flint_clz(n0);
    dcnt = flint_clz(d0);

    dbi = d0 << dcnt >> (FLINT_BITS - DIV1_NBITS);
    inv = div1_tab[dbi - (1 << (DIV1_NBITS - 1))];
    q0 = ((n0 << ncnt) >> (DIV1_NBITS + 1)) * inv >> (FLINT_BITS - 1 + ncnt - dcnt);
    r0 = n0 - q0 * d0;
    mask = -(mp_limb_t) (r0 >= d0);
    q0 -= mask;
    r0 -= d0 & mask;

    if (FLINT_UNLIKELY(r0 >= d0))
    {
        q0 = n0 / d0;
        r0 = n0 - q0 * d0;
    }

    res.d1 = q0;
    res.d0 = r0;
    return res;
}

#else
# error Unknown FLINT_MPN_HGCD2_DIV1_METHOD
#endif

/* Bit-wise two-limb division for small quotients (GMP's method 2). */
FLINT_FORCE_INLINE mp_limb_t
div2(mp_ptr rp, mp_limb_t n1, mp_limb_t n0, mp_limb_t d1, mp_limb_t d0)
{
    mp_limb_t q = 0;
    int ncnt, dcnt;

    ncnt = flint_clz(n1);
    dcnt = flint_clz(d1);
    dcnt -= ncnt;

    d1 = (d1 << dcnt) + (d0 >> 1 >> (FLINT_BITS - 1 - dcnt));
    d0 <<= dcnt;

    do
    {
        mp_limb_t mask;
        q <<= 1;
        if (FLINT_UNLIKELY(n1 == d1))
            mask = -(mp_limb_t) (n0 >= d0);
        else
            mask = -(mp_limb_t) (n1 > d1);

        q -= mask;

        sub_ddmmss(n1, n0, n1, n0, mask & d1, mask & d0);

        d0 = (d1 << (FLINT_BITS - 1)) | (d0 >> 1);
        d1 = d1 >> 1;
    }
    while (dcnt--);

    rp[0] = n0;
    rp[1] = n1;

    return q;
}

/* Reduces a, b until |a - b| (almost) fits in one limb + 1 bit. Constructs
   the matrix M. Returns 1 if we make progress, i.e. can perform at least
   one subtraction. Otherwise returns zero. */
int
_flint_mpn_hgcd2(mp_limb_t ah, mp_limb_t al, mp_limb_t bh, mp_limb_t bl,
    flint_mpn_hgcd_matrix1_struct * M)
{
    mp_limb_t u00, u01, u10, u11;

    if (ah < 2 || bh < 2)
        return 0;

    if (ah > bh || (ah == bh && al > bl))
    {
        sub_ddmmss(ah, al, ah, al, bh, bl);
        if (ah < 2)
            return 0;

        u00 = u01 = u11 = 1;
        u10 = 0;
    }
    else
    {
        sub_ddmmss(bh, bl, bh, bl, ah, al);
        if (bh < 2)
            return 0;

        u00 = u10 = u11 = 1;
        u01 = 0;
    }

    if (ah < bh)
        goto subtract_a;

    for (;;)
    {
        FLINT_ASSERT(ah >= bh);
        if (ah == bh)
            goto done;

        if (ah < (UWORD(1) << (FLINT_BITS / 2)))
        {
            ah = (ah << (FLINT_BITS / 2)) + (al >> (FLINT_BITS / 2));
            bh = (bh << (FLINT_BITS / 2)) + (bl >> (FLINT_BITS / 2));

            break;
        }

        /* Subtract a -= q b, and multiply M from the right by (1 q ; 0 1),
           affecting the second column of M. */
        FLINT_ASSERT(ah > bh);
        sub_ddmmss(ah, al, ah, al, bh, bl);

        if (ah < 2)
            goto done;

        if (ah <= bh)
        {
            /* Use q = 1 */
            u01 += u00;
            u11 += u10;
        }
        else
        {
            mp_limb_t r[2];
            mp_limb_t q = div2(r, ah, al, bh, bl);
            al = r[0]; ah = r[1];
            if (ah < 2)
            {
                /* A is too small, but q is correct. */
                u01 += q * u00;
                u11 += q * u10;
                goto done;
            }
            q++;
            u01 += q * u00;
            u11 += q * u10;
        }
subtract_a:
        FLINT_ASSERT(bh >= ah);
        if (ah == bh)
            goto done;

        if (bh < (UWORD(1) << (FLINT_BITS / 2)))
        {
            ah = (ah << (FLINT_BITS / 2)) + (al >> (FLINT_BITS / 2));
            bh = (bh << (FLINT_BITS / 2)) + (bl >> (FLINT_BITS / 2));

            goto subtract_a1;
        }

        /* Subtract b -= q a, and multiply M from the right by (1 0 ; q 1),
           affecting the first column of M. */
        sub_ddmmss(bh, bl, bh, bl, ah, al);

        if (bh < 2)
            goto done;

        if (bh <= ah)
        {
            /* Use q = 1 */
            u00 += u01;
            u10 += u11;
        }
        else
        {
            mp_limb_t r[2];
            mp_limb_t q = div2(r, bh, bl, ah, al);
            bl = r[0]; bh = r[1];
            if (bh < 2)
            {
                /* B is too small, but q is correct. */
                u00 += q * u01;
                u10 += q * u11;
                goto done;
            }
            q++;
            u00 += q * u01;
            u10 += q * u11;
        }
    }

    /* NOTE: Since we discard the least significant half limb, we don't get a
       truly maximal M (corresponding to |a - b| < 2^{FLINT_BITS + 1}). */
    /* Single precision loop */
    for (;;)
    {
        FLINT_ASSERT(ah >= bh);

        ah -= bh;
        if (ah < (UWORD(1) << (FLINT_BITS / 2 + 1)))
            break;

        if (ah <= bh)
        {
            /* Use q = 1 */
            u01 += u00;
            u11 += u10;
        }
        else
        {
            _double_limb_t rq = div1(ah, bh);
            mp_limb_t q = rq.d1;
            ah = rq.d0;

            if (ah < (UWORD(1) << (FLINT_BITS / 2 + 1)))
            {
                /* A is too small, but q is correct. */
                u01 += q * u00;
                u11 += q * u10;
                break;
            }
            q++;
            u01 += q * u00;
            u11 += q * u10;
        }
subtract_a1:
        FLINT_ASSERT(bh >= ah);

        bh -= ah;
        if (bh < (UWORD(1) << (FLINT_BITS / 2 + 1)))
            break;

        if (bh <= ah)
        {
            /* Use q = 1 */
            u00 += u01;
            u10 += u11;
        }
        else
        {
            _double_limb_t rq = div1(bh, ah);
            mp_limb_t q = rq.d1;
            bh = rq.d0;

            if (bh < (UWORD(1) << (FLINT_BITS / 2 + 1)))
            {
                /* B is too small, but q is correct. */
                u00 += q * u01;
                u10 += q * u11;
                break;
            }
            q++;
            u00 += q * u01;
            u10 += q * u11;
        }
    }

done:
    M->u[0][0] = u00; M->u[0][1] = u01;
    M->u[1][0] = u10; M->u[1][1] = u11;

    return 1;
}

/*
    The 2 x 2 single-limb matrix-vector products of the Lehmer steps.
    GMP forms each output with a mpn_mul_1 and a mpn_addmul_1 or
    mpn_submul_1 pass; here both outputs are formed in a single pass,
    which reads every input limb once and allows the first output to
    overwrite the first input (in-place use avoids the copies GMP's
    callers make). The matrix entries are less than 2^(FLINT_BITS - 1),
    so the double-limb accumulators cannot overflow.
*/

/* Sets (r;b) = (a;b) M, with M = (u00, u01; u10, u11), i.e.
   (r, b) = (u00 a + u10 b, u01 a + u11 b). r and b must have space for
   n + 1 limbs. r may be equal to a. */
mp_size_t
_flint_mpn_hgcd_mul_matrix1_vector(const flint_mpn_hgcd_matrix1_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
{
    mp_limb_t u00 = M->u[0][0], u01 = M->u[0][1];
    mp_limb_t u10 = M->u[1][0], u11 = M->u[1][1];
    mp_limb_t c0 = 0, c1 = 0;
    mp_size_t i;

    for (i = 0; i < n; i++)
    {
        mp_limb_t a = ap[i], b = bp[i], p1, p0, q1, q0;

        umul_ppmm(p1, p0, u00, a);
        umul_ppmm(q1, q0, u10, b);
        add_ssaaaa(p1, p0, p1, p0, q1, q0);
        add_ssaaaa(p1, p0, p1, p0, 0, c0);
        rp[i] = p0;
        c0 = p1;

        umul_ppmm(p1, p0, u01, a);
        umul_ppmm(q1, q0, u11, b);
        add_ssaaaa(p1, p0, p1, p0, q1, q0);
        add_ssaaaa(p1, p0, p1, p0, 0, c1);
        bp[i] = p0;
        c1 = p1;
    }

    rp[n] = c0;
    bp[n] = c1;

    n += (c0 | c1) > 0;
    return n;
}

/* Sets (r;b) = M^{-1}(a;b), with M^{-1} = (u11, -u01; -u10, u00) from the
   left, i.e. (r, b) = (u11 a - u01 b, u00 b - u10 a); the results are
   known to be nonnegative and to fit in n limbs. r may be equal to a.

   With ~x = B^n - 1 - x (the limbwise complement), u11 a - u01 b =
   u11 a + u01 ~b + u01 - u01 B^n, so each output is an unsigned sum of
   two products as in _flint_mpn_hgcd_mul_matrix1_vector, with the carry
   starting at u01 (respectively u10); since the result fits in n limbs,
   the final carry is u01 (respectively u10). */
mp_size_t
_flint_mpn_hgcd_mul_matrix1_inverse_vector(const flint_mpn_hgcd_matrix1_struct * M,
    mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
{
    mp_limb_t u00 = M->u[0][0], u01 = M->u[0][1];
    mp_limb_t u10 = M->u[1][0], u11 = M->u[1][1];
    mp_limb_t c0 = u01, c1 = u10;
    mp_size_t i;

    for (i = 0; i < n; i++)
    {
        mp_limb_t a = ap[i], b = bp[i], p1, p0, q1, q0;

        umul_ppmm(p1, p0, u11, a);
        umul_ppmm(q1, q0, u01, ~b);
        add_ssaaaa(p1, p0, p1, p0, q1, q0);
        add_ssaaaa(p1, p0, p1, p0, 0, c0);
        rp[i] = p0;
        c0 = p1;

        umul_ppmm(p1, p0, u00, b);
        umul_ppmm(q1, q0, u10, ~a);
        add_ssaaaa(p1, p0, p1, p0, q1, q0);
        add_ssaaaa(p1, p0, p1, p0, 0, c1);
        bp[i] = p0;
        c1 = p1;
    }

    FLINT_ASSERT(c0 == u01);
    FLINT_ASSERT(c1 == u10);

    n -= (rp[n - 1] | bp[n - 1]) == 0;
    return n;
}

#if !FLINT_HAVE_NATIVE_mpn_gcd_22

/* Port of GMP's generic mpn_gcd_22, used when GMP's (possibly assembly)
   version is not accessible. Both inputs must be odd. */
void
_flint_mpn_gcd_22(mp_ptr g, mp_limb_t u1, mp_limb_t u0, mp_limb_t v1, mp_limb_t v0)
{
    FLINT_ASSERT(u0 & v0 & 1);

    /* Implicit least significant bit */
    u0 = (u0 >> 1) | (u1 << (FLINT_BITS - 1));
    u1 >>= 1;

    v0 = (v0 >> 1) | (v1 << (FLINT_BITS - 1));
    v1 >>= 1;

    while (u1 || v1) /* u1 == 0 can happen at most twice per call */
    {
        mp_limb_t vgtu, t1, t0;
        sub_ddmmss(t1, t0, u1, u0, v1, v0);
        vgtu = -(t1 >> (FLINT_BITS - 1));

        if (FLINT_UNLIKELY(t0 == 0))
        {
            int c;
            if (t1 == 0)
            {
                g[1] = (u1 << 1) | (u0 >> (FLINT_BITS - 1));
                g[0] = (u0 << 1) | 1;
                return;
            }
            c = flint_ctz(t1);

            /* v1 = min (u1, v1) */
            v1 += (vgtu & t1);
            /* u0 = |u1 - v1| */
            u0 = (t1 ^ vgtu) - vgtu;
            FLINT_ASSERT(c < FLINT_BITS - 1);
            u0 >>= c + 1;
            u1 = 0;
        }
        else
        {
            int c;
            c = flint_ctz(t0);
            c++;
            /* V <-- min (U, V). */
            add_ssaaaa(v1, v0, v1, v0, vgtu & t1, vgtu & t0);
            /* U  <--  |U - V|. No carry handling needed in this conditional
               negation, since t0 != 0. */
            u0 = (t0 ^ vgtu) - vgtu;
            u1 = t1 ^ vgtu;
            if (FLINT_UNLIKELY(c == FLINT_BITS))
            {
                u0 = u1;
                u1 = 0;
            }
            else
            {
                u0 = (u0 >> c) | (u1 << (FLINT_BITS - c));
                u1 >>= c;
            }
        }
    }
    while ((v0 | u0) & (UWORD(1) << (FLINT_BITS - 1)))
    {
        /* At most two iterations */
        mp_limb_t vgtu, t0;
        int c;
        sub_ddmmss(vgtu, t0, 0, u0, 0, v0);
        if (FLINT_UNLIKELY(t0 == 0))
        {
            g[1] = u0 >> (FLINT_BITS - 1);
            g[0] = (u0 << 1) | 1;
            return;
        }

        /* v <-- min (u, v) */
        v0 += (vgtu & t0);

        /* u <-- |u - v| */
        u0 = (t0 ^ vgtu) - vgtu;

        c = flint_ctz(t0);
        u0 = (u0 >> 1) >> c;
    }

    g[0] = mpn_gcd_11((u0 << 1) + 1, (v0 << 1) + 1);
    g[1] = 0;
}

#endif

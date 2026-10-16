/*
    From GMP's mpn/generic/sbpi1_bdiv_q.c:
    Copyright 2005, 2006, 2009, 2011, 2012, 2017 Free Software Foundation, Inc.
    Contributed to the GNU project by Niels Möller and Torbjörn Granlund.

    From GMP's mpn/generic/sbpi1_bdiv_qr.c:
    Copyright 2006, 2009, 2011, 2012, 2017 Free Software Foundation, Inc.
    Contributed to the GNU project by Niels Möller and Torbjörn Granlund.

    From GMP's mpn/generic/dcpi1_bdiv_qr.c:
    Copyright 2006, 2007, 2009, 2010, 2017 Free Software Foundation, Inc.
    Contributed to the GNU project by Niels Möller and Torbjorn Granlund.

    From GMP's mpn/generic/dcpi1_bdiv_q.c:
    Copyright 2006, 2007, 2009-2011, 2017 Free Software Foundation, Inc.
    Contributed to the GNU project by Niels Möller and Torbjorn Granlund.

    From GMP's mpn/generic/divis.c:
    Copyright 2001, 2002, 2005, 2009, 2014, 2017, 2018 Free Software
    Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    Exact division by Hensel (2-adic) division: with b = 2^v B^k b' for odd
    b', a = 2^v B^k a' and q = a' / b' = a' b'^(-1) mod B^n where n is the
    number of quotient limbs. Requires that b divides a.
*/
void
_flint_mpn_divexact_hensel(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t k = 0, n;
    unsigned int v;
    mp_ptr as, bs;
    TMP_INIT;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    n = an - bn + 1;

    /* the quotient has fewer limbs when the top limb of a is below that of b */
    if (n > 1 && a[an - 1] < b[bn - 1])
    {
        q[n - 1] = 0;
        n--;
    }

    /* strip zero low limbs */
    while (b[k] == 0)
        k++;
    FLINT_ASSERT(flint_mpn_zero_p(a, k));
    a += k;
    an -= k;
    b += k;
    bn -= k;

    v = flint_ctz(b[0]);

    /* q = a' / b' mod B^n only depends on the low n limbs of a' and b'
       (one more limb of a and b is read for the shift) */
    if (an > n + 1)
        an = n + 1;
    if (bn > n + 1)
        bn = n + 1;

    if (v == 0)
    {
        flint_mpn_bdiv_q(q, a, FLINT_MIN(an, n), b, FLINT_MIN(bn, n), n);
        return;
    }

    TMP_START;
    as = TMP_ALLOC((an + bn) * sizeof(mp_limb_t));
    bs = as + an;
    mpn_rshift(as, a, an, v);
    mpn_rshift(bs, b, bn, v);
    flint_mpn_bdiv_q(q, as, FLINT_MIN(an, n), bs, FLINT_MIN(bn, n), n);
    TMP_END;
}

/*
    A bidirectional variant (low half of the quotient by Hensel division,
    high half by Euclidean division, so that two half-size problems replace
    one full-size one) was implemented and measured: it was about 1.2x
    slower than the pure Hensel division for balanced operands and 2.5x
    slower when the quotient is shorter than the divisor, since the
    Euclidean half (an inverse plus three half-size products) costs more
    than the Hensel half (an inverse plus two). Exact division is therefore
    purely 2-adic.
*/

/*
    Hensel division with a 1-limb inverse, below the Newton range.

    * Truncated divisors of at most FLINT_MPN_DIVEXACT_SMALL_BN limbs
      (_bdiv_small): the quotient limbs q_i = w b^(-1) mod B are developed
      from the bottom with the divisor and a window of the partial
      remainder in registers, the operands shifted on the fly (to make the
      divisor odd), and the product q_i b subtracted with carry chains. Only
      the low qn limbs of a and the low min(bn, qn) limbs of b matter.

    * Otherwise ports of GMP's mpn_sbpi1_bdiv_q, mpn_sbpi1_bdiv_qr and
      mpn_dcpi1_bdiv_* (Copyright 2006-2012 Free Software Foundation, Inc.,
      by Niels Moller and Torbjorn Granlund), with FLINT's multiplication.
      Following GMP, these compute Q = -N / D mod B^n with
      dinv = -1 / D mod B.
*/

/* {W, D} = ({x, W[D-1], ..., W[1]} - ({q d, D + 1} >> FLINT_BITS)
   - pend B^(D-1)), setting bw to the borrow (0 or 1), for D = 1, ..., 7 */
#define BDIV_SMALL_STEP_1(bw) \
    do { \
        mp_limb_t __h0, __l0, __r0, __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        (void) __l0; \
        sub_ddmmss(__z, __r0, 0, x, 0, __h0 + pend); \
        W[0] = __r0; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_2(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __r0, __r1, __z, __p[3]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        (void) __l0; \
        add_ssaaaa(__p[2], __p[1], \
            __h1, __h0, \
            0, __l1); \
        sub_dddmmmsss(__z, __r1, __r0, \
            0, x, W[1], \
            0, __p[2] + pend, __p[1]); \
        W[0] = __r0; W[1] = __r1; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_3(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __r0, __r1, __r2, __z, __p[4]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        (void) __l0; \
        add_sssaaaaaa(__p[3], __p[2], __p[1], \
            __h2, __h1, __h0, \
            0, __l2, __l1); \
        sub_ddddmmmmssss(__z, __r2, __r1, __r0, \
            0, x, W[2], W[1], \
            0, __p[3] + pend, __p[2], __p[1]); \
        W[0] = __r0; W[1] = __r1; W[2] = __r2; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_4(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __r0, __r1, __r2, __r3, __z, __p[5]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        (void) __l0; \
        add_ssssaaaaaaaa(__p[4], __p[3], __p[2], __p[1], \
            __h3, __h2, __h1, __h0, \
            0, __l3, __l2, __l1); \
        sub_dddddmmmmmsssss(__z, __r3, __r2, __r1, __r0, \
            0, x, W[3], W[2], W[1], \
            0, __p[4] + pend, __p[3], __p[2], __p[1]); \
        W[0] = __r0; W[1] = __r1; W[2] = __r2; W[3] = __r3; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_5(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __h4, __l4, __r0, __r1, __r2, __r3, __r4, __z, __p[6]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        umul_ppmm(__h4, __l4, q, d[4]); \
        (void) __l0; \
        add_sssssaaaaaaaaaa(__p[5], __p[4], __p[3], __p[2], __p[1], \
            __h4, __h3, __h2, __h1, __h0, \
            0, __l4, __l3, __l2, __l1); \
        sub_ddddddmmmmmmssssss(__z, __r4, __r3, __r2, __r1, __r0, \
            0, x, W[4], W[3], W[2], W[1], \
            0, __p[5] + pend, __p[4], __p[3], __p[2], __p[1]); \
        W[0] = __r0; W[1] = __r1; W[2] = __r2; W[3] = __r3; W[4] = __r4; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_6(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __h4, __l4, __h5, __l5, __r0, __r1, __r2, __r3, __r4, __r5, __z, __p[7]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        umul_ppmm(__h4, __l4, q, d[4]); \
        umul_ppmm(__h5, __l5, q, d[5]); \
        (void) __l0; \
        add_ssssssaaaaaaaaaaaa(__p[6], __p[5], __p[4], __p[3], __p[2], __p[1], \
            __h5, __h4, __h3, __h2, __h1, __h0, \
            0, __l5, __l4, __l3, __l2, __l1); \
        sub_dddddddmmmmmmmsssssss(__z, __r5, __r4, __r3, __r2, __r1, __r0, \
            0, x, W[5], W[4], W[3], W[2], W[1], \
            0, __p[6] + pend, __p[5], __p[4], __p[3], __p[2], __p[1]); \
        W[0] = __r0; W[1] = __r1; W[2] = __r2; W[3] = __r3; W[4] = __r4; W[5] = __r5; \
        (bw) = -__z; \
    } while (0)

#define BDIV_SMALL_STEP_7(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __h4, __l4, __h5, __l5, __h6, __l6, __r0, __r1, __r2, __r3, __r4, __r5, __r6, __z, __p[8]; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        umul_ppmm(__h4, __l4, q, d[4]); \
        umul_ppmm(__h5, __l5, q, d[5]); \
        umul_ppmm(__h6, __l6, q, d[6]); \
        (void) __l0; \
        add_sssssssaaaaaaaaaaaaaa(__p[7], __p[6], __p[5], __p[4], __p[3], __p[2], __p[1], \
            __h6, __h5, __h4, __h3, __h2, __h1, __h0, \
            0, __l6, __l5, __l4, __l3, __l2, __l1); \
        sub_ddddddddmmmmmmmmssssssss(__z, __r6, __r5, __r4, __r3, __r2, __r1, __r0, \
            0, x, W[6], W[5], W[4], W[3], W[2], W[1], \
            0, __p[7] + pend, __p[6], __p[5], __p[4], __p[3], __p[2], __p[1]); \
        W[0] = __r0; W[1] = __r1; W[2] = __r2; W[3] = __r3; W[4] = __r4; W[5] = __r5; W[6] = __r6; \
        (bw) = -__z; \
    } while (0)

/* q = (a / 2^v) / (b / 2^v) mod B^qn for an exact division, using the low
   D = min(bn, qn) limbs of the shifted divisor; requires an > qn (or v = 0
   and an >= qn) and 2 <= D <= FLINT_MPN_DIVEXACT_SMALL_BN */
FLINT_FORCE_INLINE void
_bdiv_small(mp_ptr qp, mp_srcptr ap, mp_size_t qn, mp_srcptr bp, mp_size_t bn,
    unsigned int v, const int D)
{
    mp_limb_t d[FLINT_MPN_DIVEXACT_SMALL_BN], W[FLINT_MPN_DIVEXACT_SMALL_BN];
    mp_limb_t binv, q, x, pend;
    unsigned int t = FLINT_BITS - 1 - v;
    mp_size_t i;
    int j;

#define SHR(lo, hi) (((lo) >> v) | (((hi) << 1) << t))

    for (j = 0; j < D; j++)
        d[j] = (j + 1 < bn) ? SHR(bp[j], bp[j + 1]) : (bp[j] >> v);

    for (j = 0; j < D; j++)
        W[j] = (v == 0) ? ap[j] : SHR(ap[j], ap[j + 1]);

    binv = n_binvert(d[0]);
    pend = 0;

    for (i = 0; ; i++)
    {
        q = W[0] * binv;
        qp[i] = q;

        if (i == qn - 1)
            break;

        /* dividend limbs at positions >= qn do not affect the quotient */
        if (i + D < qn)
            x = (v == 0) ? ap[i + D] : SHR(ap[i + D], ap[i + D + 1]);
        else
            x = 0;

        if (D == 2)
            BDIV_SMALL_STEP_2(pend);
        else if (D == 3)
            BDIV_SMALL_STEP_3(pend);
        else if (D == 4)
            BDIV_SMALL_STEP_4(pend);
        else if (D == 5)
            BDIV_SMALL_STEP_5(pend);
        else if (D == 6)
            BDIV_SMALL_STEP_6(pend);
        else
            BDIV_SMALL_STEP_7(pend);
    }

#undef SHR
}

static void
bdiv_small(mp_ptr qp, mp_srcptr ap, mp_size_t qn, mp_srcptr bp, mp_size_t bn, unsigned int v)
{
    switch (FLINT_MIN(bn, qn))
    {
        case 1:
            /* qn = 1 (bn >= 2): the low limb of a' times the inverse of
               the low limb of b' */
            FLINT_ASSERT(qn == 1 && bn >= 2);
            {
                unsigned int t = FLINT_BITS - 1 - v;
                mp_limb_t a0 = (v == 0) ? ap[0] : ((ap[0] >> v) | ((ap[1] << 1) << t));
                mp_limb_t b0 = (bp[0] >> v) | ((bp[1] << 1) << t);
                qp[0] = a0 * n_binvert(b0);
            }
            break;
        case 2: _bdiv_small(qp, ap, qn, bp, bn, v, 2); break;
        case 3: _bdiv_small(qp, ap, qn, bp, bn, v, 3); break;
        case 4: _bdiv_small(qp, ap, qn, bp, bn, v, 4); break;
        case 5: _bdiv_small(qp, ap, qn, bp, bn, v, 5); break;
        case 6: _bdiv_small(qp, ap, qn, bp, bn, v, 6); break;
        default: _bdiv_small(qp, ap, qn, bp, bn, v, 7); break;
    }
}

/* {p, n} += c with carry propagation, no carry out */
static void
_incr_u(mp_ptr p, mp_limb_t c)
{
    p[0] += c;
    if (p[0] < c)
        while (++(*++p) == 0)
            ;
}

/* Divisibility of a by b for short odd parts b' = b / 2^v (bn' <= 7
   limbs), with the Hensel remainder computed in registers: with a' =
   a / 2^v (an' limbs) and k = an' - bn' + 1 (or k = an' - bn' when the
   top limb of a' is below that of b', so that a' < b' B^k), the quotient
   q = a' / b' mod B^k is developed as in _bdiv_small, but with the full
   divisor and all limbs of a', leaving R = (a' - q b') / B^k = W - pend
   B^D with -B^D < R < B^D. b' divides a' iff its quotient (of at most k
   limbs) is q, i.e. iff R = 0. The low v bits of a must be zero. */
FLINT_FORCE_INLINE int
_divisible_small(mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn,
    unsigned int v, const int D)
{
    mp_limb_t d[FLINT_MPN_DIVEXACT_SMALL_BN], W[FLINT_MPN_DIVEXACT_SMALL_BN];
    mp_limb_t binv, q, x, pend, any;
    unsigned int t = FLINT_BITS - 1 - v;
    mp_size_t i, k;
    int j;

#define SHR(lo, hi) (((lo) >> v) | (((hi) << 1) << t))
#define A(jj) (((jj) < an) ? SHR(ap[jj], ((jj) + 1 < an) ? ap[(jj) + 1] : 0) : 0)

    for (j = 0; j < D; j++)
        d[j] = SHR(bp[j], (j + 1 < bn) ? bp[j + 1] : 0);

    /* k = an' - bn' + 1 quotient limbs for a' of an' limbs, one less when
       the top limb of a' is below that of b' (as in GMP) */
    k = an - ((ap[an - 1] >> v) == 0);
    if (k < D)
        return 0;
    k = k - D + (A(k - 1) >= d[D - 1]);
    if (k == 0)
        return 0;

    for (j = 0; j < D; j++)
        W[j] = A(j);

    binv = n_binvert(d[0]);
    pend = 0;

    for (i = 0; i < k; i++)
    {
        q = W[0] * binv;
        x = A(i + D);

        if (D == 1)
            BDIV_SMALL_STEP_1(pend);
        else if (D == 2)
            BDIV_SMALL_STEP_2(pend);
        else if (D == 3)
            BDIV_SMALL_STEP_3(pend);
        else if (D == 4)
            BDIV_SMALL_STEP_4(pend);
        else if (D == 5)
            BDIV_SMALL_STEP_5(pend);
        else if (D == 6)
            BDIV_SMALL_STEP_6(pend);
        else
            BDIV_SMALL_STEP_7(pend);
    }

    any = pend;
    for (j = 0; j < D; j++)
        any |= W[j];

#undef A
#undef SHR

    return any == 0;
}

/* requires bn >= 2, b[bn - 1] != 0, v = ctz(b[0]), the low v bits of a
   zero and a[an - 1] != 0, and b / 2^v of at most
   FLINT_MPN_DIVEXACT_SMALL_BN limbs */
int
_flint_mpn_divisible_small(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, unsigned int v)
{
    FLINT_ASSERT(bn >= 2);
    FLINT_ASSERT(a[an - 1] != 0);

    switch (bn - ((b[bn - 1] >> v) == 0))
    {
        case 1: return _divisible_small(a, an, b, bn, v, 1);
        case 2: return _divisible_small(a, an, b, bn, v, 2);
        case 3: return _divisible_small(a, an, b, bn, v, 3);
        case 4: return _divisible_small(a, an, b, bn, v, 4);
        case 5: return _divisible_small(a, an, b, bn, v, 5);
        case 6: return _divisible_small(a, an, b, bn, v, 6);
        default: return _divisible_small(a, an, b, bn, v, 7);
    }
}

/* GMP's mpn_sbpi1_bdiv_q: {qp, un} = -{up, un} / {dp, dn} mod B^un,
   destroying {up, un}; un >= dn >= 1 */
static void
bdiv_q_basecase_preinv1(mp_ptr qp, mp_ptr up, mp_size_t un, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_size_t i;
    mp_limb_t q;

    FLINT_ASSERT(dn > 0);
    FLINT_ASSERT(un >= dn);
    FLINT_ASSERT(dp[0] & 1);

    if (un > dn)
    {
        mp_limb_t cy, hi;

        for (i = un - dn - 1, cy = 0; i > 0; i--)
        {
            q = dinv * up[0];
            hi = mpn_addmul_1(up, dp, dn, q);
            *qp++ = q;
            hi += cy;
            cy = hi < cy;
            hi += up[dn];
            cy += hi < up[dn];
            up[dn] = hi;
            up++;
        }
        q = dinv * up[0];
        hi = cy + mpn_addmul_1(up, dp, dn, q);
        *qp++ = q;
        up[dn] += hi;
        up++;
    }

    for (i = dn; i > 1; i--)
    {
        q = dinv * up[0];
        mpn_addmul_1(up, dp, i, q);
        *qp++ = q;
        up++;
    }

    *qp = dinv * up[0];
}

/* GMP's mpn_sbpi1_bdiv_qr: the un - dn quotient limbs and, in the top dn
   limbs of {up, un}, the remainder; returns the carry. un > dn >= 1 */
static mp_limb_t
bdiv_qr_basecase_preinv1(mp_ptr qp, mp_ptr up, mp_size_t un, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_size_t i;
    mp_limb_t cy;

    FLINT_ASSERT(dn > 0);
    FLINT_ASSERT(un > dn);
    FLINT_ASSERT(dp[0] & 1);

    for (i = un - dn, cy = 0; i != 0; i--)
    {
        mp_limb_t q = dinv * up[0];
        mp_limb_t hi = mpn_addmul_1(up, dp, dn, q);
        *qp++ = q;
        hi += cy;
        cy = hi < cy;
        hi += up[dn];
        cy += hi < up[dn];
        up[dn] = hi;
        up++;
    }

    return cy;
}

/* {tp, xn + yn} = {x, xn} {y, yn} in either order */
static void
_mul_any(mp_ptr tp, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    if (xn >= yn)
        flint_mpn_mul(tp, x, xn, y, yn);
    else
        flint_mpn_mul(tp, y, yn, x, xn);
}

/* GMP's mpn_dcpi1_bdiv_qr_n: n quotient limbs, {np, 2n} replaced by the
   remainder in its high half; tp has n limbs */
static mp_limb_t
bdiv_qr_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t lo, hi;
    mp_limb_t cy, rh;

    lo = n >> 1;
    hi = n - lo;

    if (lo < FLINT_MPN_DC_BDIV_QR_CUTOFF)
        cy = bdiv_qr_basecase_preinv1(qp, np, 2 * lo, dp, lo, dinv);
    else
        cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, lo, dinv, tp);

    _mul_any(tp, dp + lo, hi, qp, lo);
    _incr_u(tp + lo, cy);
    rh = mpn_add(np + lo, np + lo, n + hi, tp, n);

    if (hi < FLINT_MPN_DC_BDIV_QR_CUTOFF)
        cy = bdiv_qr_basecase_preinv1(qp + lo, np + lo, 2 * hi, dp, hi, dinv);
    else
        cy = bdiv_qr_n_divconquer_preinv1(qp + lo, np + lo, dp, hi, dinv, tp);

    _mul_any(tp, qp + lo, hi, dp + hi, lo);
    _incr_u(tp + hi, cy);
    rh += mpn_add_n(np + n, np + n, tp, n);

    return rh;
}

/* GMP's mpn_dcpi1_bdiv_q_n: Q = -N / D mod B^n, destroying {np, n} */
static void
bdiv_q_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
{
    while (n >= FLINT_MPN_DC_BDIV_Q_CUTOFF)
    {
        mp_size_t lo, hi;
        mp_limb_t cy;

        lo = n >> 1;
        hi = n - lo;

        cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, lo, dinv, tp);

        flint_mpn_mullow_n(tp, qp, dp + hi, lo);
        mpn_add_n(np + hi, np + hi, tp, lo);

        if (lo < hi)
        {
            cy += mpn_addmul_1(np + lo, qp, lo, dp[lo]);
            np[n - 1] += cy;
        }

        qp += lo;
        np += lo;
        n -= lo;
    }

    bdiv_q_basecase_preinv1(qp, np, n, dp, n, dinv);
}

/* GMP's mpn_dcpi1_bdiv_q: Q = -N / D mod B^nn, destroying {np, nn};
   nn >= dn >= FLINT_MPN_DC_BDIV_Q_CUTOFF (the only use) */
static void
bdiv_q_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_size_t qn;
    mp_limb_t cy;
    mp_ptr tp;
    TMP_INIT;

    FLINT_ASSERT(dn >= FLINT_MPN_DC_BDIV_Q_CUTOFF);
    FLINT_ASSERT(nn - dn >= 0);
    FLINT_ASSERT(dp[0] & 1);

    TMP_START;
    tp = TMP_ALLOC(dn * sizeof(mp_limb_t));

    qn = nn;

    if (qn > dn)
    {
        do
            qn -= dn;
        while (qn > dn);

        if (qn < FLINT_MPN_DC_BDIV_QR_CUTOFF)
            cy = bdiv_qr_basecase_preinv1(qp, np, 2 * qn, dp, qn, dinv);
        else
            cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, qn, dinv, tp);

        if (qn != dn)
        {
            _mul_any(tp, qp, qn, dp + qn, dn - qn);
            _incr_u(tp + qn, cy);
            mpn_add(np + qn, np + qn, nn - qn, tp, dn);
            cy = 0;
        }

        np += qn;
        qp += qn;

        qn = nn - qn;
        while (qn > dn)
        {
            mpn_add_1(np + dn, np + dn, qn - dn, cy);
            cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, dn, dinv, tp);
            qp += dn;
            np += dn;
            qn -= dn;
        }

        bdiv_q_n_divconquer_preinv1(qp, np, dp, dn, dinv, tp);
    }
    else
    {
        bdiv_q_n_divconquer_preinv1(qp, np, dp, qn, dinv, tp);
    }

    TMP_END;
}

/* GMP's mpn_dcpi1_bdiv_qr: the nn - dn quotient limbs of -N / D and, in
   the top dn limbs of {np, nn}, the remainder; returns the carry.
   nn - dn and dn >= FLINT_MPN_DC_BDIV_QR_CUTOFF (the only use) */
static mp_limb_t
bdiv_qr_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_size_t qn;
    mp_limb_t rr, cy;
    mp_ptr tp;
    TMP_INIT;

    FLINT_ASSERT(dn >= FLINT_MPN_DC_BDIV_QR_CUTOFF);
    FLINT_ASSERT(nn - dn >= FLINT_MPN_DC_BDIV_QR_CUTOFF);
    FLINT_ASSERT(dp[0] & 1);

    TMP_START;
    tp = TMP_ALLOC(dn * sizeof(mp_limb_t));

    qn = nn - dn;

    if (qn > dn)
    {
        /* reduce qn mod dn, doing the typically smaller block first */
        do
            qn -= dn;
        while (qn > dn);

        if (qn < FLINT_MPN_DC_BDIV_QR_CUTOFF)
            cy = bdiv_qr_basecase_preinv1(qp, np, 2 * qn, dp, qn, dinv);
        else
            cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, qn, dinv, tp);

        rr = 0;
        if (qn != dn)
        {
            _mul_any(tp, qp, qn, dp + qn, dn - qn);
            _incr_u(tp + qn, cy);
            rr = mpn_add(np + qn, np + qn, nn - qn, tp, dn);
            cy = 0;
        }

        np += qn;
        qp += qn;

        qn = nn - dn - qn;
        do
        {
            rr += mpn_add_1(np + dn, np + dn, qn, cy);
            cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, dn, dinv, tp);
            qp += dn;
            np += dn;
            qn -= dn;
        }
        while (qn > 0);

        TMP_END;
        return rr + cy;
    }

    cy = bdiv_qr_n_divconquer_preinv1(qp, np, dp, qn, dinv, tp);

    rr = 0;
    if (qn != dn)
    {
        _mul_any(tp, qp, qn, dp + qn, dn - qn);
        _incr_u(tp + qn, cy);
        rr = mpn_add(np + qn, np + qn, nn - qn, tp, dn);
        cy = 0;
    }

    TMP_END;
    return rr + cy;
}

/* Divisibility test following GMP's mpn_divisible_p. Requires
   an >= bn >= 2, a[an - 1] != 0, b[bn - 1] != 0, v = ctz(b[0]) and the
   low v bits of a zero. With a' = a / 2^v and b' = b / 2^v (odd), padded
   so that the top limb of a' is below that of b' (whence a' < b' B^k),
   the Hensel division a' + Q b' = R B^k by k = len(a') - len(b') quotient
   limbs gives a remainder R < 2 b', and b' divides a' iff it divides R,
   i.e. iff R = b' (R = 0 is impossible as a' != 0). */
int
_flint_mpn_divisible_bdiv(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, unsigned int v)
{
    mp_ptr rp, qp, dp;
    mp_size_t k;
    mp_limb_t dinv;
    int res;
    TMP_INIT;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 2);
    FLINT_ASSERT(a[an - 1] != 0);
    FLINT_ASSERT(b[bn - 1] != 0);

    TMP_START;
    rp = TMP_ALLOC((2 * an + 2) * sizeof(mp_limb_t));
    qp = rp + an + 1;

    if (v != 0)
    {
        dp = TMP_ALLOC(bn * sizeof(mp_limb_t));
        mpn_rshift(dp, b, bn, v);
        mpn_rshift(rp, a, an, v);
        bn -= (dp[bn - 1] == 0);
        an -= (rp[an - 1] == 0);

        if (bn == 1)
        {
            res = flint_mpn_divisible_1_odd(rp, an, dp[0]);
            TMP_END;
            return res;
        }

        /* the shift can make a' shorter than b' (an = bn on entry with a
           top limb of a below 2^v): then 0 < a' < b' */
        if (an < bn)
        {
            TMP_END;
            return 0;
        }
    }
    else
    {
        dp = (mp_ptr) b;
        flint_mpn_copyi(rp, a, an);
    }

    if (rp[an - 1] >= dp[bn - 1])
    {
        rp[an] = 0;
        an++;
    }
    else if (an == bn)
    {
        TMP_END;
        return 0;
    }

    k = an - bn;
    dinv = -n_binvert(dp[0]);

    /* the returned carry is zero if R = b' (mod B^bn) */
    if (bn < FLINT_MPN_DC_BDIV_QR_CUTOFF || k < FLINT_MPN_DC_BDIV_QR_CUTOFF)
        bdiv_qr_basecase_preinv1(qp, rp, an, dp, bn, dinv);
    else
        bdiv_qr_divconquer_preinv1(qp, rp, an, dp, bn, dinv);

    res = flint_mpn_equal_p(rp + k, dp, bn);

    TMP_END;
    return res;
}

/* exact quotient q = a' / b' mod B^qn of the shifted operands (b' odd) by
   the 1-limb-inverse Hensel division; dn = min(bn, qn) >= 2 */
static void
bdiv_q_preinv1(mp_ptr qp, mp_srcptr ap, mp_size_t qn, mp_srcptr bp, mp_size_t bn, unsigned int v)
{
    mp_ptr np, dp;
    mp_size_t dn = FLINT_MIN(bn, qn);
    mp_limb_t dinv;
    TMP_INIT;

    TMP_START;
    np = TMP_ALLOC((qn + 1 + dn + 1) * sizeof(mp_limb_t));

    if (v != 0)
    {
        dp = np + qn + 1;
        /* the low qn limbs of a' and dn limbs of b' (reading one more
           limb of b when it is truncated) */
        mpn_rshift(np, ap, qn + 1, v);
        mpn_rshift(dp, bp, FLINT_MIN(bn, dn + 1), v);
    }
    else
    {
        dp = (mp_ptr) bp;
        flint_mpn_copyi(np, ap, qn);
    }

    dinv = -n_binvert(dp[0]);

    if (dn < FLINT_MPN_DC_BDIV_Q_CUTOFF)
        bdiv_q_basecase_preinv1(qp, np, qn, dp, dn, dinv);
    else
        bdiv_q_divconquer_preinv1(qp, np, qn, dp, dn, dinv);

    /* the quotient -N/D computed above, negated */
    mpn_neg(qp, qp, qn);

    TMP_END;
}

/* Exact division for bn >= 1 (b divides a); the top limb of a may be
   zero. Only the low limbs of the operands are used. */
void
_flint_mpn_divexact(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
{
    mp_size_t qn, dn;
    unsigned int v;

    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(b[bn - 1] != 0);

    /* strip zero low limbs */
    while (b[0] == 0)
    {
        FLINT_ASSERT(a[0] == 0);
        a++;
        an--;
        b++;
        bn--;
    }

    qn = an - bn + 1;

    if (bn == 1)
    {
        mpn_divexact_1(q, a, an, b[0]);
        return;
    }

    /* the top quotient limb vanishes when the top limb of a is below that
       of b */
    if (qn > 1 && a[an - 1] < b[bn - 1])
    {
        q[qn - 1] = 0;
        qn--;
    }

    v = flint_ctz(b[0]);
    dn = FLINT_MIN(bn, qn);

    if (dn <= FLINT_MPN_DIVEXACT_SMALL_BN)
        bdiv_small(q, a, qn, b, bn, v);
    else if (dn >= FLINT_MPN_DIVEXACT_NEWTON_CUTOFF
            || (bn >= FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF && qn >= 4 * bn))
        _flint_mpn_divexact_hensel(q, a, an, b, bn);
    else
        bdiv_q_preinv1(q, a, qn, b, bn, v);
}

/*
    From GMP's mpn/generic/sbpi1_div_qr.c, sbpi1_divappr_q.c, sbpi1_div_q.c:
    Copyright 2007, 2009 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund.

    From GMP's mpn/generic/dcpi1_div_qr.c:
    Copyright 2006, 2007, 2009 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund.

    From GMP's mpn/generic/dcpi1_divappr_q.c:
    Copyright 2006, 2007, 2009, 2010 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund.

    From GMP's mpn/generic/dcpi1_div_q.c:
    Copyright 2006, 2007, 2009, 2010 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund and Marco Bodrato.

    From GMP's mpn/generic/div_q.c:
    Copyright 2009, 2010, 2015, 2018 Free Software Foundation, Inc.
    Contributed to the GNU project by Torbjorn Granlund.

    From GMP's mpn/generic/tdiv_qr.c:
    Copyright 1997, 2000-2002, 2005, 2009, 2015 Free Software Foundation, Inc.

    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"
#include "mp_real.h"

/*
    Euclidean and approximate division below the Newton range, and the
    dispatcher _flint_mpn_tdiv_qr. The functions call each other through
    static versions (the exported _flint_mpn_* names are wrappers at the end
    of the file), avoiding the call overhead of exported symbols in shared
    libraries.

    * Divisors of at most FLINT_MPN_DIV_SMALL_BN limbs (_div_small): one
      function per divisor length, with the (normalized) divisor and the
      partial remainder kept in registers and the dividend normalized on the
      fly. The algorithm is the schoolbook division of GMP's
      mpn_sbpi1_div_qr (a 3/2 division by the top two divisor limbs per
      quotient limb, followed by a subtraction of the product of the
      quotient limb and the remaining divisor limbs), without function calls
      or memory traffic. For one-limb divisors, a chain of hardware
      divisions is used for short dividends if these are fast
      (FLINT_PREINVERT_LIMB_USE_NATIVE); otherwise GMP's assembly
      mpn_divrem_1.

    * Division with a normalized divisor and a precomputed 3/2 inverse,
      ported from GMP (see the copyright notices above): the schoolbook
      functions mpn_sbpi1_*, the divide and conquer
      functions mpn_dcpi1_* (using FLINT's multiplication) and the top-level
      logic of mpn_tdiv_qr and mpn_div_q (using FLINT's schoolbook, divide
      and conquer and Newton division and FLINT's multiplication).
*/

/* static versions of the exported functions */
static void tdiv_qr_small(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn);
static mp_limb_t divrem_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
static mp_limb_t divapprox_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
static mp_limb_t div_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv);
static mp_limb_t divrem_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp);
static mp_limb_t divapprox_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp);
static mp_limb_t divrem_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp);
static mp_limb_t divapprox_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp);
static void tdiv_qr_divconquer(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn);
static void tdiv_q_divconquer(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn, int exact);
FLINT_FORCE_INLINE void divapprox(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn);

/* {R, BN} = (R B + x) - (B - 1) d, which is < d, when the top two limbs of
   R and d agree (the case excluded by the 3/2 division) */
static void
_div_small_step_special(mp_limb_t * R, mp_limb_t x, const mp_limb_t * d, int bn)
{
    mp_limb_t t[FLINT_MPN_DIV_SMALL_BN + 1], c, c1, u;
    int i;

    /* t = R B + x - d B + d, computed mod B^(bn + 1) */
    t[0] = x;
    for (i = 0; i < bn; i++)
        t[i + 1] = R[i];

    c = 0;
    for (i = 0; i < bn; i++)
    {
        u = t[i + 1] - d[i];
        c1 = t[i + 1] < d[i];
        c1 += u < c;
        t[i + 1] = u - c;
        c = c1;
    }

    c = 0;
    for (i = 0; i < bn; i++)
    {
        u = t[i] + d[i];
        c1 = u < d[i];
        u += c;
        c1 += u < c;
        t[i] = u;
        c = c1;
    }

    for (i = 0; i < bn; i++)
        R[i] = t[i];
}

/* {n1, n0, R[BN-4], ..., R[0], x} -= q {d, BN - 2}, setting bw to the
   borrow (0 or -1), for BN = 3, ..., 7 */
#define DIV_SMALL_STEP_3(bw) \
    do { \
        mp_limb_t __h0, __l0, __p[2], __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        __p[0] = __l0; \
        __p[1] = __h0; \
        sub_ddddmmmmssss(__z, n1, n0, x, \
            0, n1, n0, x, \
            0, 0, __p[1], __p[0]); \
        (bw) = __z; \
    } while (0)

#define DIV_SMALL_STEP_4(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __p[3], __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        add_sssaaaaaa(__p[2], __p[1], __p[0], \
            __h1, __h0, 0, \
            0, __l1, __l0); \
        sub_dddddmmmmmsssss(__z, n1, n0, R[0], x, \
            0, n1, n0, R[0], x, \
            0, 0, __p[2], __p[1], __p[0]); \
        (bw) = __z; \
    } while (0)

#define DIV_SMALL_STEP_5(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __p[4], __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        add_ssssaaaaaaaa(__p[3], __p[2], __p[1], __p[0], \
            __h2, __h1, __h0, 0, \
            0, __l2, __l1, __l0); \
        sub_ddddddmmmmmmssssss(__z, n1, n0, R[1], R[0], x, \
            0, n1, n0, R[1], R[0], x, \
            0, 0, __p[3], __p[2], __p[1], __p[0]); \
        (bw) = __z; \
    } while (0)

#define DIV_SMALL_STEP_6(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __p[5], __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        add_sssssaaaaaaaaaa(__p[4], __p[3], __p[2], __p[1], __p[0], \
            __h3, __h2, __h1, __h0, 0, \
            0, __l3, __l2, __l1, __l0); \
        sub_dddddddmmmmmmmsssssss(__z, n1, n0, R[2], R[1], R[0], x, \
            0, n1, n0, R[2], R[1], R[0], x, \
            0, 0, __p[4], __p[3], __p[2], __p[1], __p[0]); \
        (bw) = __z; \
    } while (0)

#define DIV_SMALL_STEP_7(bw) \
    do { \
        mp_limb_t __h0, __l0, __h1, __l1, __h2, __l2, __h3, __l3, __h4, __l4, __p[6], __z; \
        umul_ppmm(__h0, __l0, q, d[0]); \
        umul_ppmm(__h1, __l1, q, d[1]); \
        umul_ppmm(__h2, __l2, q, d[2]); \
        umul_ppmm(__h3, __l3, q, d[3]); \
        umul_ppmm(__h4, __l4, q, d[4]); \
        add_ssssssaaaaaaaaaaaa(__p[5], __p[4], __p[3], __p[2], __p[1], __p[0], \
            __h4, __h3, __h2, __h1, __h0, 0, \
            0, __l4, __l3, __l2, __l1, __l0); \
        sub_ddddddddmmmmmmmmssssssss(__z, n1, n0, R[3], R[2], R[1], R[0], x, \
            0, n1, n0, R[3], R[2], R[1], R[0], x, \
            0, 0, __p[5], __p[4], __p[3], __p[2], __p[1], __p[0]); \
        (bw) = __z; \
    } while (0)

/* HW: use hardware 2/1 divisions (FLINT_MPN_UDIV_QR_3BY2_HW) instead of
   the 3/2 inverse, saving its computation for short dividends */
FLINT_FORCE_INLINE void
_div_small(mp_ptr qp, mp_ptr rp, mp_srcptr ap, const mp_size_t an, mp_srcptr bp, const int BN, const int HW)
{
    mp_limb_t d[FLINT_MPN_DIV_SMALL_BN], R[FLINT_MPN_DIV_SMALL_BN], top;
    unsigned int s, t;
    mp_size_t i, start;
    int j;

    FLINT_ASSERT(an >= BN);
    FLINT_ASSERT(bp[BN - 1] != 0);

    /* x << s | y >> (FLINT_BITS - s), also for s = 0 */
    s = flint_clz(bp[BN - 1]);
    t = FLINT_BITS - 1 - s;

#define SHL(x, y) (((x) << s) | (((y) >> 1) >> t))
#define NLIMB(k) (((k) == 0) ? (ap[0] << s) : SHL(ap[k], ap[(k) - 1]))

    for (j = BN - 1; j > 0; j--)
        d[j] = SHL(bp[j], bp[j - 1]);
    d[0] = bp[0] << s;

    /* The shifted dividend has an + 1 limbs, the top one being smaller than
       d[BN - 1], so that the top BN limbs are smaller than d. If the top
       limb is zero (always the case for a normalized divisor), the top
       quotient limb is 0 or 1 and follows from a comparison instead, as in
       GMP; this saves one step. */
    top = (ap[an - 1] >> 1) >> t;
    if (top != 0)
    {
        R[BN - 1] = top;
        for (j = BN - 2; j >= 0; j--)
            R[j] = NLIMB(an - BN + 1 + j);
        start = an - BN;
    }
    else
    {
        mp_limb_t T[FLINT_MPN_DIV_SMALL_BN], bw, u;

        for (j = BN - 1; j >= 0; j--)
            R[j] = NLIMB(an - BN + j);

        /* R >= d ? R - d : R */
        bw = 0;
        for (j = 0; j < BN; j++)
        {
            u = R[j] - d[j];
            T[j] = u - bw;
            bw = (R[j] < d[j]) | (u < bw);
        }
        if (!bw)
            for (j = 0; j < BN; j++)
                R[j] = T[j];

        qp[an - BN] = !bw;
        start = an - BN - 1;
    }

    if (BN == 1)
    {
        mp_limb_t v, r = R[0];

        FLINT_MPN_INVERT_LIMB(v, d[0]);
        for (i = start; i >= 0; i--)
            FLINT_MPN_UDIV_QR_2BY1(qp[i], r, r, NLIMB(i), d[0], v);
        R[0] = r;
    }
    else if (BN == 2)
    {
        mp_limb_t v, r1 = R[1], r0 = R[0];

        v = flint_mpn_preinv1(d[1], d[0]);
        for (i = start; i >= 0; i--)
            FLINT_MPN_UDIV_QR_3BY2(qp[i], r1, r0, r1, r0, NLIMB(i), d[1], d[0], v);
        R[0] = r0;
        R[1] = r1;
    }
    else
    {
        mp_limb_t v = 0, q, n1, n0, x, bw, u, c, c1;

        if (!HW)
            v = flint_mpn_preinv1(d[BN - 1], d[BN - 2]);

        for (i = start; i >= 0; i--)
        {
            x = NLIMB(i);

            if (FLINT_UNLIKELY(R[BN - 1] == d[BN - 1] && R[BN - 2] == d[BN - 2]))
            {
                _div_small_step_special(R, x, d, BN);
                qp[i] = UWORD_MAX;
                continue;
            }

            if (HW)
                FLINT_MPN_UDIV_QR_3BY2_HW(q, n1, n0, R[BN - 1], R[BN - 2], R[BN - 3], d[BN - 1], d[BN - 2]);
            else
                FLINT_MPN_UDIV_QR_3BY2(q, n1, n0, R[BN - 1], R[BN - 2], R[BN - 3], d[BN - 1], d[BN - 2], v);

            if (BN == 3)
                DIV_SMALL_STEP_3(bw);
            else if (BN == 4)
                DIV_SMALL_STEP_4(bw);
            else if (BN == 5)
                DIV_SMALL_STEP_5(bw);
            else if (BN == 6)
                DIV_SMALL_STEP_6(bw);
            else
                DIV_SMALL_STEP_7(bw);

            /* the new remainder is {n1, n0, R[BN-4], ..., R[0], x} */
            for (j = BN - 3; j >= 1; j--)
                R[j] = R[j - 1];
            R[0] = x;
            R[BN - 2] = n0;
            R[BN - 1] = n1;

            if (FLINT_UNLIKELY(bw != 0))
            {
                c = 0;
                for (j = 0; j < BN; j++)
                {
                    u = R[j] + d[j];
                    c1 = u < d[j];
                    u += c;
                    c1 += u < c;
                    R[j] = u;
                    c = c1;
                }
                q--;
            }

            qp[i] = q;
        }
    }

    if (rp != NULL)
    {
        for (j = 0; j < BN - 1; j++)
            rp[j] = (R[j] >> s) | ((R[j + 1] << 1) << t);
        rp[BN - 1] = R[BN - 1] >> s;
    }

#undef SHL
#undef NLIMB
}

#if FLINT_PREINVERT_LIMB_USE_NATIVE
/* one-limb divisor with hardware division, no normalization needed */
static mp_limb_t
_div_small_1_hw(mp_ptr qp, mp_srcptr ap, mp_size_t an, mp_limb_t b)
{
    mp_limb_t r;
    mp_size_t i;

    r = ap[an - 1];
    if (r < b)
    {
        qp[an - 1] = 0;
    }
    else
    {
        qp[an - 1] = r / b;
        r = r % b;
    }

    for (i = an - 2; i >= 0; i--)
        udiv_qrnnd(qp[i], r, r, ap[i], b);

    return r;
}
#endif

#if FLINT_PREINVERT_LIMB_USE_NATIVE
/* two-limb divisor by hardware divisions, without an inverse (see
   FLINT_MPN_UDIV_QR_3BY2_HW); the dividend is normalized on the fly */
static void
_div_small_2_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp)
{
    mp_limb_t d1, d0, r1, r0, n;
    unsigned int s, t;
    mp_size_t i;

    s = flint_clz(bp[1]);
    t = FLINT_BITS - 1 - s;

    /* x << s for the high part hi and the next lower limb lo, without an
       undefined shift for s = 0 */
#define SHL(hi, lo) (((hi) << s) | (((lo) >> 1) >> t))
    d1 = SHL(bp[1], bp[0]);
    d0 = bp[0] << s;
    r1 = (ap[an - 1] >> 1) >> t;
    r0 = SHL(ap[an - 1], (an >= 2) ? ap[an - 2] : 0);

    for (i = an - 2; i >= 0; i--)
    {
        n = SHL(ap[i], (i >= 1) ? ap[i - 1] : 0);
        FLINT_MPN_UDIV_QR_3BY2_HW(qp[i], r1, r0, r1, r0, n, d1, d0);
    }
#undef SHL

    if (rp != NULL)
    {
        rp[0] = (r0 >> s) | ((r1 << 1) << t);
        rp[1] = r1 >> s;
    }
}
#endif

#if FLINT_HAVE_NATIVE_mpn_divrem_2
mp_limb_t __gmpn_divrem_2(mp_ptr, mp_size_t, mp_ptr, mp_size_t, mp_srcptr);

/* two-limb divisor by GMP's assembly mpn_divrem_2, which needs a
   normalized divisor and overwrites the dividend: the dividend is shifted
   into a temporary */
static void
_div_2_gmp(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp)
{
    mp_limb_t d[2];
    mp_ptr t;
    unsigned int s;
    TMP_INIT;

    TMP_START;
    t = TMP_ALLOC((an + 1) * sizeof(mp_limb_t));
    s = flint_clz(bp[1]);

    if (s != 0)
    {
        d[1] = (bp[1] << s) | (bp[0] >> (FLINT_BITS - s));
        d[0] = bp[0] << s;
        /* t[an] < 2^s <= d[1], so the quotient has an - 1 limbs */
        t[an] = mpn_lshift(t, ap, an, s);
        __gmpn_divrem_2(qp, 0, t, an + 1, d);
        if (rp != NULL)
        {
            rp[0] = (t[0] >> s) | (t[1] << (FLINT_BITS - s));
            rp[1] = t[1] >> s;
        }
    }
    else
    {
        flint_mpn_copyi(t, ap, an);
        qp[an - 2] = __gmpn_divrem_2(qp, 0, t, an, bp);
        if (rp != NULL)
        {
            rp[0] = t[0];
            rp[1] = t[1];
        }
    }

    TMP_END;
}
#endif

static void _div_small_2(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 2, 0); }
static void _div_small_3(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 3, 0); }
static void _div_small_4(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 4, 0); }
static void _div_small_5(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 5, 0); }
static void _div_small_6(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 6, 0); }
static void _div_small_7(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 7, 0); }
#if FLINT_PREINVERT_LIMB_USE_NATIVE
static void _div_small_3_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 3, 1); }
static void _div_small_4_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 4, 1); }
static void _div_small_5_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 5, 1); }
static void _div_small_6_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 6, 1); }
static void _div_small_7_hw(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp) { _div_small(qp, rp, ap, an, bp, 7, 1); }
#endif

static void
tdiv_qr_small(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    FLINT_ASSERT(an >= bn);
    FLINT_ASSERT(bn >= 1 && bn <= FLINT_MPN_DIV_SMALL_BN);
    FLINT_ASSERT(bp[bn - 1] != 0);

    switch (bn)
    {
        case 1:
        {
            mp_limb_t r;
            if (an == 1)
            {
                qp[0] = ap[0] / bp[0];
                r = ap[0] % bp[0];
            }
            else
            {
#if FLINT_PREINVERT_LIMB_USE_NATIVE
                /* the hardware division chain is latency bound; GMP's
                   pipelined mpn_divrem_1 wins for long dividends */
                if (FLINT_MPN_DIVREM_1_USE_HW(an, bp[0]))
                    r = _div_small_1_hw(qp, ap, an, bp[0]);
                else
#endif
                    r = mpn_divrem_1(qp, 0, ap, an, bp[0]);
            }
            if (rp != NULL)
                rp[0] = r;
            break;
        }
        case 2:
#if FLINT_PREINVERT_LIMB_USE_NATIVE
            if (an < FLINT_MPN_DIV_2_HW_CUTOFF)
            {
                _div_small_2_hw(qp, rp, ap, an, bp);
                break;
            }
#endif
#if FLINT_HAVE_NATIVE_mpn_divrem_2
            if (an >= FLINT_MPN_DIV_2_GMP_CUTOFF)
            {
                _div_2_gmp(qp, rp, ap, an, bp);
                break;
            }
#endif
            _div_small_2(qp, rp, ap, an, bp);
            break;
        default:
#if FLINT_PREINVERT_LIMB_USE_NATIVE
            /* short quotients: hardware divisions save computing the 3/2
               inverse, but each step is slower */
            if (an - bn + 1 < FLINT_MPN_DIV_SMALL_HW_QN_CUTOFF)
            {
                switch (bn)
                {
                    case 3: _div_small_3_hw(qp, rp, ap, an, bp); break;
                    case 4: _div_small_4_hw(qp, rp, ap, an, bp); break;
                    case 5: _div_small_5_hw(qp, rp, ap, an, bp); break;
                    case 6: _div_small_6_hw(qp, rp, ap, an, bp); break;
                    default: _div_small_7_hw(qp, rp, ap, an, bp); break;
                }
                break;
            }
#endif
            switch (bn)
            {
                case 3: _div_small_3(qp, rp, ap, an, bp); break;
                case 4: _div_small_4(qp, rp, ap, an, bp); break;
                case 5: _div_small_5(qp, rp, ap, an, bp); break;
                case 6: _div_small_6(qp, rp, ap, an, bp); break;
                default: _div_small_7(qp, rp, ap, an, bp); break;
            }
            break;
    }
}

#define FLINT_UDIV_QR_3BY2 FLINT_MPN_UDIV_QR_3BY2

/* Schoolbook division, ported from GMP's mpn_sbpi1_div_qr:
   {np, nn} = {qp, nn - dn} {dp, dn} + {np, dn}, returning the high
   quotient limb. Requires dn > 2, nn >= dn, dp normalized. */
static mp_limb_t
divrem_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn,
    mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_limb_t qh, n1, n0, d1, d0, cy, cy1, q;
    mp_size_t i;

    FLINT_ASSERT(dn > 2);
    FLINT_ASSERT(nn >= dn);
    FLINT_ASSERT(dp[dn - 1] >> (FLINT_BITS - 1));

    np += nn;

    qh = mpn_cmp(np - dn, dp, dn) >= 0;
    if (qh)
        mpn_sub_n(np - dn, np - dn, dp, dn);

    qp += nn - dn;

    dn -= 2;
    d1 = dp[dn + 1];
    d0 = dp[dn + 0];

    np -= 2;
    n1 = np[1];

    for (i = nn - (dn + 2); i > 0; i--)
    {
        np--;
        if (FLINT_UNLIKELY(n1 == d1) && np[1] == d0)
        {
            q = UWORD_MAX;
            mpn_submul_1(np - dn, dp, dn + 2, q);
            n1 = np[1];
        }
        else
        {
            FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

            cy = mpn_submul_1(np - dn, dp, dn, q);

            cy1 = n0 < cy;
            n0 = n0 - cy;
            cy = n1 < cy1;
            n1 = n1 - cy1;
            np[0] = n0;

            if (FLINT_UNLIKELY(cy != 0))
            {
                n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
                q--;
            }
        }

        *--qp = q;
    }

    np[1] = n1;
    return qh;
}

/* Approximate schoolbook division, ported from GMP's mpn_sbpi1_divappr_q:
   {qp, nn - dn} and the returned high limb form a quotient that is correct
   or one too large; {np, nn} is destroyed. Requires dn > 2, nn >= dn, dp
   normalized. */
static mp_limb_t
divapprox_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn,
    mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_limb_t qh, n1, n0, d1, d0, cy, cy1, q, flag;
    mp_size_t qn, i;

    FLINT_ASSERT(dn > 2);
    FLINT_ASSERT(nn >= dn);
    FLINT_ASSERT(dp[dn - 1] >> (FLINT_BITS - 1));

    np += nn;

    qn = nn - dn;
    if (qn + 1 < dn)
    {
        dp += dn - (qn + 1);
        dn = qn + 1;
    }

    qh = mpn_cmp(np - dn, dp, dn) >= 0;
    if (qh)
        mpn_sub_n(np - dn, np - dn, dp, dn);

    qp += qn;

    dn -= 2;
    d1 = dp[dn + 1];
    d0 = dp[dn + 0];

    np -= 2;
    n1 = np[1];

    for (i = qn - (dn + 2); i >= 0; i--)
    {
        np--;
        if (FLINT_UNLIKELY(n1 == d1) && np[1] == d0)
        {
            q = UWORD_MAX;
            mpn_submul_1(np - dn, dp, dn + 2, q);
            n1 = np[1];
        }
        else
        {
            FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

            cy = mpn_submul_1(np - dn, dp, dn, q);

            cy1 = n0 < cy;
            n0 = n0 - cy;
            cy = n1 < cy1;
            n1 -= cy1;
            np[0] = n0;

            if (FLINT_UNLIKELY(cy != 0))
            {
                n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
                q--;
            }
        }

        *--qp = q;
    }

    flag = ~UWORD(0);

    if (dn >= 0)
    {
        for (i = dn; i > 0; i--)
        {
            np--;
            if (FLINT_UNLIKELY(n1 >= (d1 & flag)))
            {
                q = UWORD_MAX;
                cy = mpn_submul_1(np - dn, dp, dn + 2, q);

                if (FLINT_UNLIKELY(n1 != cy))
                {
                    if (n1 < (cy & flag))
                    {
                        q--;
                        mpn_add_n(np - dn, np - dn, dp, dn + 2);
                    }
                    else
                        flag = 0;
                }
                n1 = np[1];
            }
            else
            {
                FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

                cy = mpn_submul_1(np - dn, dp, dn, q);

                cy1 = n0 < cy;
                n0 = n0 - cy;
                cy = n1 < cy1;
                n1 -= cy1;
                np[0] = n0;

                if (FLINT_UNLIKELY(cy != 0))
                {
                    n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
                    q--;
                }
            }

            *--qp = q;

            /* truncate the operands */
            dn--;
            dp++;
        }

        np--;
        if (FLINT_UNLIKELY(n1 >= (d1 & flag)))
        {
            q = UWORD_MAX;
            cy = mpn_submul_1(np, dp, 2, q);

            if (FLINT_UNLIKELY(n1 != cy))
            {
                if (n1 < (cy & flag))
                {
                    q--;
                    add_ssaaaa(np[1], np[0], np[1], np[0], dp[1], dp[0]);
                }
                else
                    flag = 0;
            }
            n1 = np[1];
        }
        else
        {
            FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

            np[1] = n1;
            np[0] = n0;
        }

        *--qp = q;
    }

    FLINT_ASSERT(np[1] == n1);

    return qh;
}

/* Exact schoolbook quotient, ported from GMP's mpn_sbpi1_div_q:
   {qp, nn - dn} and the returned high limb form floor({np, nn} / {dp, dn});
   {np, nn} is destroyed. Like divapprox_basecase_preinv1, the last dn - 2
   quotient limbs are developed with a truncated divisor; the result is then
   corrected (rarely needed) using the ignored parts. Requires dn > 2,
   nn >= dn, dp normalized. */
static mp_limb_t
div_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn,
    mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    mp_limb_t qh, n1, n0, d1, d0, cy, cy1, q, flag;
    mp_size_t qn, i;
    mp_size_t dn_orig = dn;
    mp_srcptr dp_orig = dp;
    mp_ptr np_orig = np;

    FLINT_ASSERT(dn > 2);
    FLINT_ASSERT(nn >= dn);
    FLINT_ASSERT(dp[dn - 1] >> (FLINT_BITS - 1));

    np += nn;

    qn = nn - dn;
    if (qn + 1 < dn)
    {
        dp += dn - (qn + 1);
        dn = qn + 1;
    }

    qh = mpn_cmp(np - dn, dp, dn) >= 0;
    if (qh)
        mpn_sub_n(np - dn, np - dn, dp, dn);

    qp += qn;

    dn -= 2;
    d1 = dp[dn + 1];
    d0 = dp[dn + 0];

    np -= 2;
    n1 = np[1];

    for (i = qn - (dn + 2); i >= 0; i--)
    {
        np--;
        if (FLINT_UNLIKELY(n1 == d1) && np[1] == d0)
        {
            q = UWORD_MAX;
            mpn_submul_1(np - dn, dp, dn + 2, q);
            n1 = np[1];
        }
        else
        {
            FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

            cy = mpn_submul_1(np - dn, dp, dn, q);

            cy1 = n0 < cy;
            n0 = n0 - cy;
            cy = n1 < cy1;
            n1 -= cy1;
            np[0] = n0;

            if (FLINT_UNLIKELY(cy != 0))
            {
                n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
                q--;
            }
        }

        *--qp = q;
    }

    flag = ~UWORD(0);

    if (dn >= 0)
    {
        for (i = dn; i > 0; i--)
        {
            np--;
            if (FLINT_UNLIKELY(n1 >= (d1 & flag)))
            {
                q = UWORD_MAX;
                cy = mpn_submul_1(np - dn, dp, dn + 2, q);

                if (FLINT_UNLIKELY(n1 != cy))
                {
                    if (n1 < (cy & flag))
                    {
                        q--;
                        mpn_add_n(np - dn, np - dn, dp, dn + 2);
                    }
                    else
                        flag = 0;
                }
                n1 = np[1];
            }
            else
            {
                FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

                cy = mpn_submul_1(np - dn, dp, dn, q);

                cy1 = n0 < cy;
                n0 = n0 - cy;
                cy = n1 < cy1;
                n1 -= cy1;
                np[0] = n0;

                if (FLINT_UNLIKELY(cy != 0))
                {
                    n1 += d1 + mpn_add_n(np - dn, np - dn, dp, dn + 1);
                    q--;
                }
            }

            *--qp = q;

            /* truncate the operands */
            dn--;
            dp++;
        }

        np--;
        if (FLINT_UNLIKELY(n1 >= (d1 & flag)))
        {
            q = UWORD_MAX;
            cy = mpn_submul_1(np, dp, 2, q);

            if (FLINT_UNLIKELY(n1 != cy))
            {
                if (n1 < (cy & flag))
                {
                    q--;
                    add_ssaaaa(np[1], np[0], np[1], np[0], dp[1], dp[0]);
                }
                else
                    flag = 0;
            }
            n1 = np[1];
        }
        else
        {
            FLINT_UDIV_QR_3BY2(q, n1, n0, n1, np[1], np[0], d1, d0, dinv);

            np[0] = n0;
            np[1] = n1;
        }

        *--qp = q;
    }

    FLINT_ASSERT(np[1] == n1);
    np += 2;

    dn = dn_orig;
    if (FLINT_UNLIKELY(n1 < (dn & flag)))
    {
        mp_limb_t x;

        /* the quotient may be too large if the remainder is small:
           recompute with the ignored operand parts until the remainder
           spills */
        x = n1;
        if (dn > 2)
        {
            /* compensate for the triangularization */
            mp_limb_t y;

            dp = dp_orig;
            if (qn + 1 < dn)
            {
                dp += dn - (qn + 1);
                dn = qn + 1;
            }

            y = np[-2];

            for (i = dn - 3; i >= 0; i--)
            {
                q = qp[i];
                cy = mpn_submul_1(np - (dn - i), dp, dn - i - 2, q);
                if (y < cy)
                {
                    if (x == 0)
                    {
                        cy = mpn_sub_1(qp, qp, qn, 1);
                        FLINT_ASSERT(cy == 0);
                        return qh - cy;
                    }
                    x--;
                }
                y -= cy;
            }
            np[-2] = y;
        }

        dn = dn_orig;
        if (qn + 1 < dn)
        {
            /* compensate for the ignored dividend and divisor tails */
            dp = dp_orig;
            np = np_orig;

            if (qh != 0)
            {
                cy = mpn_sub_n(np + qn, np + qn, dp, dn - (qn + 1));
                if (cy != 0)
                {
                    if (x == 0)
                    {
                        if (qn != 0)
                            cy = mpn_sub_1(qp, qp, qn, 1);
                        return qh - cy;
                    }
                    x--;
                }
            }

            if (qn == 0)
                return qh;

            for (i = dn - qn - 2; i >= 0; i--)
            {
                cy = mpn_submul_1(np + i, qp, qn, dp[i]);
                cy = mpn_sub_1(np + qn + i, np + qn + i, dn - qn - i - 1, cy);
                if (cy != 0)
                {
                    if (x == 0)
                    {
                        mpn_sub_1(qp, qp, qn, 1);
                        return qh;
                    }
                    x--;
                }
            }
        }
    }

    return qh;
}

/* division by a normalized two-limb divisor, in place */
static mp_limb_t
divrem_2_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_limb_t d1, mp_limb_t d0, mp_limb_t dinv)
{
    mp_limb_t qh, r1, r0;
    mp_size_t i;

    r1 = np[nn - 1];
    r0 = np[nn - 2];
    qh = (r1 > d1) || (r1 == d1 && r0 >= d0);
    if (qh)
        sub_ddmmss(r1, r0, r1, r0, d1, d0);
    for (i = nn - 3; i >= 0; i--)
        FLINT_UDIV_QR_3BY2(qp[i], r1, r0, r1, r0, np[i], d1, d0, dinv);
    np[1] = r1;
    np[0] = r0;
    return qh;
}

/* Division by a normalized two-limb divisor {dp, 2}, in place, as
   divrem_2_preinv1 but computing the inverse itself. Without fast
   limb inversion, GMP's assembly mpn_divrem_2 is used when available. */
#if FLINT_HAVE_NATIVE_mpn_divrem_2 && !FLINT_PREINVERT_LIMB_USE_NATIVE
mp_limb_t __gmpn_divrem_2(mp_ptr, mp_size_t, mp_ptr, mp_size_t, mp_srcptr);
# define DIV_QR_2_NORM(qp, np, nn, dp) __gmpn_divrem_2((qp), 0, (np), (nn), (dp))
#else
# define DIV_QR_2_NORM(qp, np, nn, dp) \
    divrem_2_preinv1((qp), (np), (nn), (dp)[1], (dp)[0], \
        flint_mpn_preinv1((dp)[1], (dp)[0]))
#endif

/* {tp, hi + lo} = {x, hi} {y, lo} for hi = lo or lo + 1 */
static inline void
mul_hl(mp_ptr tp, mp_srcptr x, mp_size_t hi, mp_srcptr y, mp_size_t lo)
{
    flint_mpn_mul_n(tp, x, y, lo);
    if (hi != lo)
        tp[2 * lo] = mpn_addmul_1(tp + lo, y, lo, x[lo]);
}

/* {tp, xn + yn} = {x, xn} {y, yn}, either order */
static inline void
mul_any(mp_ptr tp, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
{
    if (xn >= yn)
        flint_mpn_mul(tp, x, xn, y, yn);
    else
        flint_mpn_mul(tp, y, yn, x, xn);
}

/* {np, 2n} = ({qp, n} + qh B^n) {dp, n} + {np, n}, returning qh; port of
   GMP's mpn_dcpi1_div_qr_n. tp: n limbs */
static mp_limb_t
divrem_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n,
    mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t lo, hi;
    mp_limb_t cy, qh, ql;

    /* the schoolbook base cases need at least 3 limbs, hence cutoffs >= 6 */
    FLINT_ASSERT(n >= 6);

    lo = n >> 1;
    hi = n - lo;

    if (hi < FLINT_MPN_DIV_DC_CUTOFF)
        qh = divrem_basecase_preinv1(qp + lo, np + 2 * lo, 2 * hi, dp + lo, hi, dinv);
    else
        qh = divrem_n_divconquer_preinv1(qp + lo, np + 2 * lo, dp + lo, hi, dinv, tp);

    mul_hl(tp, qp + lo, hi, dp, lo);

    cy = mpn_sub_n(np + lo, np + lo, tp, n);
    if (qh != 0)
        cy += mpn_sub_n(np + n, np + n, dp, lo);

    while (cy != 0)
    {
        qh -= mpn_sub_1(qp + lo, qp + lo, hi, 1);
        cy -= mpn_add_n(np + lo, np + lo, dp, n);
    }

    if (lo < FLINT_MPN_DIV_DC_CUTOFF)
        ql = divrem_basecase_preinv1(qp, np + hi, 2 * lo, dp + hi, lo, dinv);
    else
        ql = divrem_n_divconquer_preinv1(qp, np + hi, dp + hi, lo, dinv, tp);

    mul_hl(tp, dp, hi, qp, lo);

    cy = mpn_sub_n(np, np, tp, n);
    if (ql != 0)
        cy += mpn_sub_n(np + lo, np + lo, dp, hi);

    while (cy != 0)
    {
        mpn_sub_1(qp, qp, lo, 1);
        cy -= mpn_add_n(np, np, dp, n);
    }

    return qh;
}

/* as divrem_n_divconquer_preinv1, but {qp, n} + qh B^n is only an approximate
   quotient and {np, 2n} is destroyed; port of GMP's mpn_dcpi1_divappr_q_n.
   The low half is divided by the top half of the divisor, so each level of
   the recursion can make the quotient a few units too large (never too
   small); the callers use it with a guard limb */
static mp_limb_t
divapprox_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n,
    mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t lo, hi;
    mp_limb_t cy, qh, ql;

    /* the schoolbook base cases need at least 3 limbs, hence cutoffs >= 6 */
    FLINT_ASSERT(n >= 6);

    lo = n >> 1;
    hi = n - lo;

    if (hi < FLINT_MPN_DIV_DC_CUTOFF)
        qh = divrem_basecase_preinv1(qp + lo, np + 2 * lo, 2 * hi, dp + lo, hi, dinv);
    else
        qh = divrem_n_divconquer_preinv1(qp + lo, np + 2 * lo, dp + lo, hi, dinv, tp);

    mul_hl(tp, qp + lo, hi, dp, lo);

    cy = mpn_sub_n(np + lo, np + lo, tp, n);
    if (qh != 0)
        cy += mpn_sub_n(np + n, np + n, dp, lo);

    while (cy != 0)
    {
        qh -= mpn_sub_1(qp + lo, qp + lo, hi, 1);
        cy -= mpn_add_n(np + lo, np + lo, dp, n);
    }

    if (lo < FLINT_MPN_DIVAPPR_DC_CUTOFF)
        ql = divapprox_basecase_preinv1(qp, np + hi, 2 * lo, dp + hi, lo, dinv);
    else
        ql = divapprox_n_divconquer_preinv1(qp, np + hi, dp + hi, lo, dinv, tp);

    if (FLINT_UNLIKELY(ql != 0))
        flint_mpn_store(qp, lo, UWORD_MAX);

    return qh;
}

/* first quotient block of qn <= dn limbs, from the top dn + qn limbs of the
   partial remainder ending at np + qn (np points qn limbs below the top),
   dp pointing at the end of the divisor; in place, returns the high limb */
static mp_limb_t
div_qr_block(mp_ptr qp, mp_ptr np, mp_size_t qn, mp_srcptr dp,
    mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
{
    mp_limb_t qh, cy;

    if (qn < FLINT_MPN_DIV_DC_CUTOFF)
        return divrem_basecase_preinv1(qp, np - dn, dn + qn, dp - dn, dn, dinv);

    /* divide the top 2qn limbs by the top qn limbs of the divisor, then
       correct with the rest of the divisor */
    qh = divrem_n_divconquer_preinv1(qp, np - qn, dp - qn, qn, dinv, tp);

    if (qn != dn)
    {
        mul_any(tp, qp, qn, dp - dn, dn - qn);
        cy = mpn_sub_n(np - dn, np - dn, tp, dn);
        if (qh != 0)
            cy += mpn_sub_n(np - dn + qn, np - dn + qn, dp - dn, dn - qn);
        while (cy != 0)
        {
            qh -= mpn_sub_1(qp, qp, qn, 1);
            cy -= mpn_add_n(np - dn, np - dn, dp - dn, dn);
        }
    }

    return qh;
}

/*
    {np, nn} = ({qp, nn - dn} + qh B^(nn - dn)) {dp, dn} + {np, dn}, returning
    qh, for normalized {dp, dn} with dn >= 2, nn >= dn and dinv =
    flint_mpn_preinv1(dp[dn - 1], dp[dn - 2]). Port of GMP's
    mpn_dcpi1_div_qr (schoolbook below the cutoffs). tp needs dn limbs.
*/
static mp_limb_t
divrem_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp,
    mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t qn = nn - dn;
    mp_limb_t qh;

    FLINT_ASSERT(dn >= 2);
    FLINT_ASSERT(nn >= dn);
    FLINT_ASSERT(dp[dn - 1] >> (FLINT_BITS - 1));

    if (dn == 2)
        return divrem_2_preinv1(qp, np, nn, dp[1], dp[0], dinv);

    if (dn < FLINT_MPN_DIV_DC_CUTOFF || qn < FLINT_MPN_DIV_DC_CUTOFF)
        return divrem_basecase_preinv1(qp, np, nn, dp, dn, dinv);

    qp += qn;
    np += nn;
    dp += dn;

    if (qn > dn)
    {
        /* reduce qn mod dn, doing the smaller block first */
        do
            qn -= dn;
        while (qn > dn);

        qp -= qn;
        np -= qn;
        qh = div_qr_block(qp, np, qn, dp, dn, dinv, tp);

        qn = nn - dn - qn;
        do
        {
            qp -= dn;
            np -= dn;
            divrem_n_divconquer_preinv1(qp, np - dn, dp - dn, dn, dinv, tp);
            qn -= dn;
        }
        while (qn > 0);
    }
    else
    {
        qp -= qn;
        np -= qn;
        qh = div_qr_block(qp, np, qn, dp, dn, dinv, tp);
    }

    return qh;
}

/*
    Approximate quotient {qp, nn - dn} + qh B^(nn - dn) (returning qh) of
    {np, nn} by {dp, dn}, which is correct or one too large; the input
    {np, nn} is destroyed. Same requirements as divrem_preinv1.
    Port of GMP's mpn_dcpi1_divappr_q (schoolbook below the cutoffs), except
    that when nn = 2 dn - 1 the numerator is padded with a zero limb
    instead of reading the limb below it.
*/
static mp_limb_t
divapprox_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp,
    mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t qn = nn - dn;
    mp_limb_t qh, qsave;

    FLINT_ASSERT(dn >= 2);
    FLINT_ASSERT(nn >= dn);
    FLINT_ASSERT(dp[dn - 1] >> (FLINT_BITS - 1));

    if (dn == 2)
        return divrem_2_preinv1(qp, np, nn, dp[1], dp[0], dinv);

    /* for long quotients, all but the last block are divisions with
       remainder, for which divide and conquer pays off earlier */
    if (qn >= dn ? dn < FLINT_MPN_DIV_DC_CUTOFF
                 : (dn < FLINT_MPN_DIVAPPR_DC_CUTOFF || qn < FLINT_MPN_DIVAPPR_DC_CUTOFF))
        return divapprox_basecase_preinv1(qp, np, nn, dp, dn, dinv);

    if (qn >= dn)
    {
        mp_ptr np0 = np;

        qp += qn;
        np += nn;
        dp += dn;

        /* pretend an extra quotient limb is needed (the guard limb) */
        qn++;
        do
            qn -= dn;
        while (qn > dn);

        qp -= qn;
        np -= qn;
        qh = div_qr_block(qp, np, qn, dp, dn, dinv, tp);

        qn = nn - dn - qn + 1;
        while (qn > dn)
        {
            qp -= dn;
            np -= dn;
            divrem_n_divconquer_preinv1(qp, np - dn, dp - dn, dn, dinv, tp);
            qn -= dn;
        }

        /* dn - 1 = qn quotient limbs are left: develop them plus a guard
           limb, which is then discarded */
        qn--;
        qp -= qn;
        np -= dn;
        FLINT_ASSERT(np - dn == np0 - 1);
        qsave = qp[qn];
        /* the low limb of the 2dn-limb numerator lies one limb below np0;
           it only affects the guard limb, so use a padded copy */
        {
            mp_ptr T;
            TMP_INIT;
            TMP_START;
            T = TMP_ALLOC((2 * dn + qn + 1) * sizeof(mp_limb_t));
            T[0] = 0;
            flint_mpn_copyi(T + 1, np0, 2 * dn - 1);
            if (dn < FLINT_MPN_DIVAPPR_DC_CUTOFF)
                divapprox_basecase_preinv1(T + 2 * dn, T, 2 * dn, dp - dn, dn, dinv);
            else
                divapprox_n_divconquer_preinv1(T + 2 * dn, T, dp - dn, dn, dinv, tp);
            flint_mpn_copyi(qp, T + 2 * dn + 1, qn);
            TMP_END;
        }
        qp[qn] = qsave;
    }
    else
    {
        mp_ptr q2p, T;
        TMP_INIT;
        TMP_START;

        /* use the top qn + 1 limbs of the divisor and the top 2 (qn + 1)
           limbs of the numerator, padded with a zero limb if dn = qn + 1 */
        q2p = TMP_ALLOC((qn + 1) * sizeof(mp_limb_t));
        if (dn >= qn + 2)
        {
            qh = divapprox_n_divconquer_preinv1(q2p, np + nn - 2 * (qn + 1),
                dp + dn - (qn + 1), qn + 1, dinv, tp);
        }
        else
        {
            T = TMP_ALLOC(2 * (qn + 1) * sizeof(mp_limb_t));
            T[0] = 0;
            flint_mpn_copyi(T + 1, np, nn);
            qh = divapprox_n_divconquer_preinv1(q2p, T, dp, qn + 1, dinv, tp);
        }
        flint_mpn_copyi(qp, q2p + 1, qn);

        TMP_END;
    }

    return qh;
}

/*
    Euclidean division by the preinverted functions, a port of GMP's
    mpn_tdiv_qr for divisors of at least two limbs (with FLINT's schoolbook,
    divide and conquer and Newton division and FLINT's multiplication).
    If the quotient is at least about as long as the divisor, the divisor
    and the dividend are shifted to normalize the divisor, and the remainder
    is shifted back. Otherwise (tdiv_qr_short) the top 2 qn limbs
    of the dividend are divided by the top qn limbs of the divisor, and the
    quotient (at most 2 too large) is corrected using the product of the
    quotient and the ignored divisor limbs.
*/

static void
tdiv_qr_full(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    mp_ptr np, dp, tp;
    mp_limb_t dinv, qh;
    unsigned int cnt;
    mp_size_t nn;
    TMP_INIT;

    TMP_START;
    np = TMP_ALLOC((An + 1 + 2 * Bn) * sizeof(mp_limb_t));
    tp = np + An + 1;
    cnt = flint_clz(B[Bn - 1]);

    if (cnt != 0)
    {
        dp = tp + Bn;
        mpn_lshift(dp, B, Bn, cnt);
        np[An] = mpn_lshift(np, A, An, cnt);
        nn = An + 1;
    }
    else
    {
        dp = (mp_ptr) B;
        flint_mpn_copyi(np, A, An);
        nn = An;
    }

    dinv = flint_mpn_preinv1(dp[Bn - 1], dp[Bn - 2]);
    qh = divrem_preinv1(Q, np, nn, dp, Bn, dinv, tp);

    if (cnt == 0)
        Q[An - Bn] = qh;
    else
    {
        FLINT_ASSERT(qh == 0);
    }

    if (cnt != 0)
        mpn_rshift(R, np, Bn, cnt);
    else
        flint_mpn_copyi(R, np, Bn);

    TMP_END;
}

/* Requires An + adjust < 2 Bn, where adjust = (A[An - 1] >= B[Bn - 1]). */
static void
tdiv_qr_short(mp_ptr qp, mp_ptr rp, mp_srcptr np, mp_size_t nn,
    mp_srcptr dp, mp_size_t dn, int adjust)
{
    mp_size_t qn, in, rn;
    mp_ptr n2p, d2p, tp, q2p, r2p, rr;
    mp_limb_t cy, quotient_too_large;
    unsigned int cnt;
    TMP_INIT;

    qn = nn - dn;
    qp[qn] = 0;
    qn += adjust;

    if (qn == 0)
    {
        flint_mpn_copyi(rp, np, dn);
        return;
    }

    /* number of (at least partially) ignored limbs of the operands */
    in = dn - qn;
    FLINT_ASSERT(in >= 1);

    TMP_START;
    d2p = TMP_ALLOC((qn + (2 * qn + 1) + 2 * (qn + 1) + dn) * sizeof(mp_limb_t));
    n2p = d2p + qn;
    q2p = n2p + 2 * qn + 1;
    r2p = q2p + qn + 1;
    tp = r2p + qn + 1;

    /* normalize the top qn limbs of the divisor and shift the top 2 qn
       limbs of the dividend by the same amount */
    cnt = flint_clz(dp[dn - 1]);
    if (cnt != 0)
    {
        mpn_lshift(d2p, dp + in, qn, cnt);
        d2p[0] |= dp[in - 1] >> (FLINT_BITS - cnt);
        cy = mpn_lshift(n2p, np + nn - 2 * qn, 2 * qn, cnt);
        if (adjust)
        {
            n2p[2 * qn] = cy;
            n2p++;
        }
        else
        {
            n2p[0] |= np[nn - 2 * qn - 1] >> (FLINT_BITS - cnt);
        }
    }
    else
    {
        d2p = (mp_ptr) dp + in;
        flint_mpn_copyi(n2p, np + nn - 2 * qn, 2 * qn);
        if (adjust)
        {
            n2p[2 * qn] = 0;
            n2p++;
        }
    }

    /* approximate quotient from the extracted operands, whose quotient
       fits in qn limbs; the remainder replaces the low qn limbs of n2p */
    if (qn == 1)
    {
        mp_limb_t q0, r0;
        udiv_qrnnd(q0, r0, n2p[1], n2p[0], d2p[0]);
        n2p[0] = r0;
        qp[0] = q0;
    }
    else if (qn == 2)
    {
        mp_limb_t qh;
        qh = DIV_QR_2_NORM(qp, n2p, 4, d2p);
        FLINT_ASSERT(qh == 0);
        (void) qh;
    }
    else if (qn < FLINT_MPN_TDIV_QR_NEWTON_CUTOFF)
    {
        /* in place, with the normalized divisor */
        mp_limb_t qh;
        qh = divrem_preinv1(qp, n2p, 2 * qn, d2p, qn,
            flint_mpn_preinv1(d2p[qn - 1], d2p[qn - 2]), q2p);
        FLINT_ASSERT(qh == 0);
        (void) qh;
    }
    else
    {
        _flint_mpn_tdiv_qr(q2p, r2p, n2p, 2 * qn, d2p, qn);
        FLINT_ASSERT(q2p[qn] == 0);
        flint_mpn_copyi(qp, q2p, qn);
        flint_mpn_copyi(n2p, r2p, qn);
    }

    rn = qn;

    /* If the product of the most significant quotient limb and the first
       ignored divisor limb exceeds the top limb of the partial remainder,
       the quotient is too large. This catches most cases where the quotient
       is too large, and all cases where it is 2 too large. */
    {
        mp_limb_t dl, x, h, dummy;

        dl = (in >= 2) ? dp[in - 2] : 0;
        x = (dp[in - 1] << cnt) | ((dl >> 1) >> (FLINT_BITS - 1 - cnt));
        umul_ppmm(h, dummy, x, qp[qn - 1]);
        (void) dummy;

        if (n2p[qn - 1] < h)
        {
            mpn_sub_1(qp, qp, qn, 1);
            cy = mpn_add_n(n2p, n2p, d2p, qn);
            if (cy)
            {
                /* the partial remainder is safely large */
                n2p[qn] = cy;
                rn++;
            }
        }
    }

    quotient_too_large = 0;

    if (cnt != 0)
    {
        mp_limb_t cy1, cy2;

        /* append the partially used dividend limb to the partial remainder */
        cy1 = mpn_lshift(n2p, n2p, rn, FLINT_BITS - cnt);
        n2p[0] |= np[in - 1] & (UWORD_MAX >> cnt);

        /* update the partial remainder with the partially used divisor limb */
        cy2 = mpn_submul_1(n2p, qp, qn, dp[in - 1] & (UWORD_MAX >> cnt));
        if (qn != rn)
        {
            FLINT_ASSERT(n2p[qn] >= cy2);
            n2p[qn] -= cy2;
        }
        else
        {
            n2p[qn] = cy1 - cy2;
            quotient_too_large = (cy1 < cy2);
            rn++;
        }
        in--;
    }

    /* The partial remainder {n2p, rn} is now unshifted; subtract the
       product of the quotient and the ignored divisor limbs. The result
       has dn limbs, but the intermediate {rr, in + rn} may have dn + 1
       limbs (only possible without shift), in which case a temporary
       buffer is used. */
    rr = (in + rn > dn) ? TMP_ALLOC((dn + 1) * sizeof(mp_limb_t)) : rp;

    if (in == 0)
    {
        FLINT_ASSERT(rn == dn);
        flint_mpn_copyi(rr, n2p, rn);
    }
    else
    {
        if (in < qn)
            flint_mpn_mul(tp, qp, qn, dp, in);
        else
            flint_mpn_mul(tp, dp, in, qp, qn);

        cy = mpn_sub(n2p, n2p, rn, tp + in, qn);
        flint_mpn_copyi(rr + in, n2p, rn);
        quotient_too_large |= cy;
        cy = mpn_sub_n(rr, np, tp, in);
        cy = mpn_sub_1(rr + in, rr + in, rn, cy);
        quotient_too_large |= cy;
    }

    if (quotient_too_large)
    {
        mpn_sub_1(qp, qp, qn, 1);
        mpn_add_n(rr, rr, dp, dn);
    }

    if (rr != rp)
        flint_mpn_copyi(rp, rr, dn);

    TMP_END;
}

static void
tdiv_qr_divconquer(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    int adjust;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 2);
    FLINT_ASSERT(B[Bn - 1] != 0);

    if (R == NULL)
    {
        tdiv_q_divconquer(Q, A, An, B, Bn, 1);
        return;
    }

    /* conservative test for the quotient size */
    adjust = A[An - 1] >= B[Bn - 1];

    if (An + adjust >= 2 * Bn)
        tdiv_qr_full(Q, R, A, An, B, Bn);
    else
        tdiv_qr_short(Q, R, A, An, B, Bn, adjust);
}

/*
    Quotient only with the divide and conquer approximate division, as
    GMP's mpn_dcpi1_div_q: an approximate quotient Q' of
    A B^1 (one guard limb; B here is the limb radix) by the divisor is
    floor(A B^1 / b) or one more. Writing floor(A B^1 / b) = Q B^1 + f,
    if the guard limb of Q' is nonzero then Q = Q' >> FLINT_BITS exactly;
    otherwise Q' >> FLINT_BITS is Q or Q + 1, and the sign of the remainder
    A - Q b (with |A - Q b| < b) decides, which only needs the low Bn + 1
    limbs of the product. Without exact, Q or Q + 1 is returned.
*/
static void
div_divconquer_preinv1(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn, int exact)
{
    mp_ptr np, dp, tp, qp;
    mp_limb_t dinv, qh;
    unsigned int cnt;
    mp_size_t nn, n = An - Bn + 1;
    TMP_INIT;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 2);
    FLINT_ASSERT(B[Bn - 1] != 0);

    TMP_START;
    np = TMP_ALLOC((An + 2 + 2 * Bn + n + 2) * sizeof(mp_limb_t));
    tp = np + An + 2;
    qp = tp + 2 * Bn;
    cnt = flint_clz(B[Bn - 1]);

    np[0] = 0;
    if (cnt != 0)
    {
        dp = tp + Bn;
        mpn_lshift(dp, B, Bn, cnt);
        np[An + 1] = mpn_lshift(np + 1, A, An, cnt);
        nn = An + 2;
    }
    else
    {
        dp = (mp_ptr) B;
        flint_mpn_copyi(np + 1, A, An);
        nn = An + 1;
    }

    dinv = flint_mpn_preinv1(dp[Bn - 1], dp[Bn - 2]);
    qh = divapprox_preinv1(qp, np, nn, dp, Bn, dinv, tp);

    /* {qp, n + 1} is the approximate quotient including the guard limb */
    if (cnt == 0)
        qp[n] = qh;
    else if (FLINT_UNLIKELY(qh != 0))
    {
        /* Q' = B^(n + 1): floor(A B / b) = B^(n + 1) - 1 */
        flint_mpn_store(Q, n, UWORD_MAX);
        TMP_END;
        return;
    }

    flint_mpn_copyi(Q, qp + 1, n);

    if (exact && qp[0] == 0)
    {
        /* Q or Q + 1: the sign of A - Q b follows from the low Bn + 1 limbs
           (np is free) */
        mp_ptr t = np;
        /* An > Bn: this is only used for quotients of at least the divide
           and conquer cutoffs */
        FLINT_ASSERT(An > Bn);
        flint_mpn_mulmid(t, B, Bn, Q, n, 0, Bn + 1);
        mpn_sub_n(t, A, t, Bn + 1);
        if ((slong) t[Bn] < 0)
            mpn_sub_1(Q, Q, n, 1);
    }

    TMP_END;
}


/*
    Short (Mulders-style) approximate division. For {X, n} <= {D, n} with D
    normalized, sets {Q, n} + qh B^n (qh the return value) to an
    approximation of X B^n / D, where the low half of the numerator is zero:
    a numerator limb below B^n changes the quotient by less than 2 units,
    so that only the top n limbs are carried. The top hi = ceil(n/2)
    quotient limbs come from an exact division by the top hi limbs of D;
    subtracting their product with the rest of D is then only needed above
    B^n, i.e. for the high half of a hi x hi product (flint_mpn_mulhigh_n),
    and the low limbs follow recursively from the resulting hi limbs and
    the top lo limbs of D. Compared with the divide and conquer approximate
    division (GMP's dcpi1_divappr_q), this replaces a full product by a
    high product at each level.

    Each level adds at most a few units to the error (the truncated
    numerator, the high product, the truncated divisor and the fractional
    parts dropped in the corrections), so that the total error is a small
    multiple of log(n) units (at most 8 in tests with random operands); the
    callers allow E = B^(1/2). dinv = flint_mpn_preinv1(D[n-1], D[n-2]).
    tp: 4 n + 16 limbs.
*/
static mp_limb_t
divapprox_short_rec(mp_ptr Q, mp_srcptr X, mp_srcptr D, mp_size_t n,
    mp_limb_t dinv, mp_ptr tp)
{
    mp_size_t lo, hi;
    mp_limb_t qh, ql, cy;
    mp_ptr N, H, Dl;

    if (n < FLINT_MAX(FLINT_MPN_DIVAPPROX_SHORT_CUTOFF / 2, 6))
    {
        /* basecase on X B^n */
        N = tp;
        flint_mpn_zero(N, n);
        flint_mpn_copyi(N + n, X, n);
        return divapprox_basecase_preinv1(Q, N, 2 * n, D, n, dinv);
    }

    lo = n / 2;
    hi = n - lo;

    /* the top 2 hi limbs of X B^n divided exactly by the top hi limbs of
       D, giving the high quotient limbs and the remainder R in {N, hi} */
    N = tp;
    if (2 * hi == n)
    {
        flint_mpn_copyi(N, X, n);
    }
    else
    {
        N[0] = 0;
        flint_mpn_copyi(N + 1, X, n);
    }
    qh = divrem_preinv1(Q + lo, N, 2 * hi, D + lo, hi, dinv, tp + 2 * hi);

    /* H = high hi limbs of Q_hi D_lo', where D_lo' = D_lo if hi = lo and
       D_lo B otherwise, so that Z = R - H - qh D_lo' is the part of the
       remainder above B^n (resp. B^(n - 1)) */
    H = tp + 2 * hi;
    if (hi == lo)
    {
        Dl = (mp_ptr) D;
    }
    else
    {
        Dl = tp + 3 * hi;
        Dl[0] = 0;
        flint_mpn_copyi(Dl + 1, D, lo);
    }
    flint_mpn_mulhigh_n(H, Q + lo, Dl, hi);

    cy = mpn_sub_n(N, N, H, hi);
    if (qh)
        cy += mpn_sub_n(N, N, Dl, hi);

    /* bring Z into [0, D_hi): each unit of Q_hi is worth D_hi (plus a
       fraction below one unit). Z < D_hi already holds: without a borrow,
       Z <= R < D_hi, and otherwise D_hi is added to a negative value */
    while (cy != 0)
    {
        qh -= mpn_sub_1(Q + lo, Q + lo, hi, 1);
        cy -= mpn_add_n(N, N, D + lo, hi);
    }
    FLINT_ASSERT(mpn_cmp(N, D + lo, hi) < 0);

    /* the low quotient limbs: Z (resp. Z / B, lo limbs) over the top lo
       limbs of D */
    ql = divapprox_short_rec(Q, N + (hi - lo), D + hi, lo, dinv, tp + hi);
    if (ql)
        flint_mpn_store(Q, lo, UWORD_MAX);

    return qh;
}

/* error bound of divapprox_short, in units of the guard limb */
#define DIVAPPROX_SHORT_ERR (UWORD(1) << (FLINT_BITS / 2))

/* {D, m} = the top m limbs of b 2^s (s = clz of the top limb), zero padded
   below if Bn < m */
static void
_top_limbs_normalized(mp_ptr D, mp_size_t m, mp_srcptr B, mp_size_t Bn, unsigned int s)
{
    mp_size_t take = FLINT_MIN(Bn, m);

    flint_mpn_zero(D, m - take);
    if (s != 0)
    {
        mpn_lshift(D + m - take, B + Bn - take, take, s);
        if (take < Bn)
            D[m - take] |= B[Bn - take - 1] >> (FLINT_BITS - s);
    }
    else
    {
        flint_mpn_copyi(D + m - take, B + Bn - take, take);
    }
}

/* Approximation Qg of floor(A B^(f + 1) / b) with qn + 1 limbs
   (qn = An + f - Bn + 1, f >= 0) and error at most DIVAPPROX_SHORT_ERR, from
   the top qn + 1 limbs of A B^f (the f zero limbs are not formed) and of b;
   returns 1 if the approximation is B^(qn + 1) or more (then
   floor(A B^f / b) = B^qn - 1), in which case Qg is undefined. */
static int
divapprox_short(mp_ptr Qg, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn,
    mp_size_t f)
{
    mp_size_t qn = An + f - Bn + 1, m = qn + 1, j, k;
    mp_ptr D, X, tp;
    mp_limb_t dinv, qh, hi, lo;
    unsigned int s;
    TMP_INIT;

    TMP_START;
    D = TMP_ALLOC((6 * m + 24) * sizeof(mp_limb_t));
    X = D + m;
    tp = X + m;

    s = flint_clz(B[Bn - 1]);
    _top_limbs_normalized(D, m, B, Bn, s);

    /* X = floor(A B^f 2^s / B^(Bn - 1)), which has at most m limbs, so
       that X B^m / D ~= A B^(f + 1) 2^s / (b 2^s); limb k of A B^f is
       A[k - f] for f <= k < An + f and zero otherwise */
#define VA(kk) (((kk) >= f && (kk) < An + f) ? A[(kk) - f] : 0)
    for (j = 0; j < m; j++)
    {
        k = Bn - 1 + j;
        hi = VA(k);
        if (s != 0)
        {
            lo = (k >= 1) ? VA(k - 1) : 0;
            hi = (hi << s) | (lo >> (FLINT_BITS - s));
        }
        X[j] = hi;
    }
#undef VA

    dinv = flint_mpn_preinv1(D[m - 1], D[m - 2]);
    qh = divapprox_short_rec(Qg, X, D, m, dinv, tp);

    TMP_END;
    return qh != 0;
}

/* sign of A - Q b for a candidate quotient {Q, qn} within one of
   floor(A / b): nonzero if negative; A must have more than Bn limbs */
static int
_div_q_check_negative(mp_srcptr A, mp_srcptr B, mp_size_t Bn,
    mp_srcptr Q, mp_size_t qn)
{
    mp_ptr t;
    int neg;
    TMP_INIT;

    /* |A - Q b| < 2 b: the low Bn + 1 limbs decide */
    TMP_START;
    t = TMP_ALLOC((Bn + 1) * sizeof(mp_limb_t));
    flint_mpn_mulmid(t, B, Bn, Q, qn, 0, Bn + 1);
    mpn_sub_n(t, A, t, Bn + 1);
    neg = ((slong) t[Bn] < 0);
    TMP_END;
    return neg;
}

/* floor(A B^f / b) (exact, only for f = 0) or floor(A B^f / b) or
   floor(A B^f / b) + 1 (!exact) by the short division, for qn <= Bn + 2;
   not inlined, to keep the stack frame of the caller small for short
   quotients */
FLINT_STATIC_NOINLINE void
tdiv_q_short(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn,
    mp_size_t f, int exact)
{
    mp_size_t qn = An + f - Bn + 1;
    mp_limb_t g;
    mp_ptr Qg;
    TMP_INIT;

    TMP_START;
    Qg = TMP_ALLOC((qn + 1) * sizeof(mp_limb_t));

    FLINT_ASSERT(f == 0 || !exact);
    /* An > Bn for the checks of the exact quotient: the short division is
       only used for quotients of more than a few limbs */
    FLINT_ASSERT(!exact || An > Bn);

    if (divapprox_short(Qg, A, An, B, Bn, f))
    {
        flint_mpn_store(Q, qn, UWORD_MAX);
        TMP_END;
        return;
    }

    /* with q = {Qg + 1, qn} and g = Qg[0]: floor(A / b) is q if
       E <= g < B - E, q or q - 1 if g < E, and q or q + 1 otherwise */
    flint_mpn_copyi(Q, Qg + 1, qn);
    g = Qg[0];

    if (g >= -DIVAPPROX_SHORT_ERR)
    {
        if (mpn_add_1(Q, Q, qn, 1))
        {
            /* floor(A / b) < B^qn */
            flint_mpn_store(Q, qn, UWORD_MAX);
        }
        else if (exact && _div_q_check_negative(A, B, Bn, Q, qn))
        {
            mpn_sub_1(Q, Q, qn, 1);
        }
    }
    else if (exact && g < DIVAPPROX_SHORT_ERR)
    {
        if (_div_q_check_negative(A, B, Bn, Q, qn))
            mpn_sub_1(Q, Q, qn, 1);
    }

    TMP_END;
}

/*
    Quotient only, a port of GMP's mpn_div_q. If the quotient is at least
    about as long as the divisor, schoolbook division with a truncated
    divisor in the last steps (div_basecase_preinv1) or the divide and
    conquer approximate division with a guard limb (div_divconquer_preinv1).
    Otherwise, the top 2 qn + 1 limbs of the dividend are divided
    approximately by the top qn + 1 limbs of the divisor, giving a quotient
    with a guard limb; unless the guard limb is small, its high part is the
    quotient, and otherwise the quotient is checked by a multiplication.
    Without exact, the checks are skipped and the quotient is floor(A/B)
    or floor(A/B) + 1.
*/
static void
tdiv_q_divconquer(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn, int exact)
{
    mp_ptr new_np, new_dp, tp, scratch;
    mp_limb_t cy, qh, dinv;
    mp_size_t new_nn, qn;
    unsigned int cnt;
    TMP_INIT;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 2);
    FLINT_ASSERT(B[Bn - 1] != 0);

    qn = An - Bn + 1;
    cnt = flint_clz(B[Bn - 1]);

    /* short division for quotients up to about as long as the divisor */
    if (qn + 1 >= FLINT_MPN_DIVAPPROX_SHORT_CUTOFF && qn <= Bn + 2)
    {
        tdiv_q_short(Q, A, An, B, Bn, 0, exact);
        return;
    }

    if (qn + 5 >= Bn)   /* GMP's FUDGE = 5 */
    {
        /* divide and conquer: for balanced operands from the tuned cutoff
           (the schoolbook division with a truncated divisor saves most
           there), for long quotients as for the division with remainder */
        if ((Bn >= FLINT_MPN_TDIV_Q_DC_CUTOFF && An - Bn >= FLINT_MPN_TDIV_Q_DC_CUTOFF)
            || (Bn >= FLINT_MPN_DIV_DC_CUTOFF && qn >= 3 * Bn)
            || (Bn >= 2 * FLINT_MPN_DIV_DC_CUTOFF && 2 * qn >= 3 * Bn))
        {
            div_divconquer_preinv1(Q, A, An, B, Bn, exact);
            return;
        }

        TMP_START;
        new_np = TMP_ALLOC((An + 1 + Bn) * sizeof(mp_limb_t));

        if (cnt != 0)
        {
            new_dp = new_np + An + 1;
            cy = mpn_lshift(new_np, A, An, cnt);
            new_np[An] = cy;
            new_nn = An + (cy != 0);
            mpn_lshift(new_dp, B, Bn, cnt);
        }
        else
        {
            cy = 0;
            new_dp = (mp_ptr) B;
            flint_mpn_copyi(new_np, A, An);
            new_nn = An;
        }

        if (Bn == 2)
        {
            qh = DIV_QR_2_NORM(Q, new_np, new_nn, new_dp);
        }
        else
        {
            dinv = flint_mpn_preinv1(new_dp[Bn - 1], new_dp[Bn - 2]);
            /* the exact version ends with a correction phase costing
               O(Bn^2) whenever the remainder is small, which an
               approximate quotient does not need */
            if (exact)
                qh = div_basecase_preinv1(Q, new_np, new_nn, new_dp, Bn, dinv);
            else
                qh = divapprox_basecase_preinv1(Q, new_np, new_nn, new_dp, Bn, dinv);
        }

        if (cy == 0)
            Q[qn - 1] = qh;
        else
        {
            FLINT_ASSERT(qh == 0);
        }

        TMP_END;
        return;
    }

    /* short quotient: fixed-size buffers for small qn */
    {
        mp_limb_t tp_s[6 * 16 + 8];   /* 6 qn + 7 limbs, including the dc scratch */
        mp_ptr rp;
        int small = (qn <= 16);

        TMP_START;
        new_nn = 2 * qn + 1;
        if (small)
            tp = tp_s;
        else
            tp = TMP_ALLOC(((qn + 2) + (new_nn + 1) + (qn + 1) + 2 * (qn + 1)) * sizeof(mp_limb_t));
        new_np = tp + qn + 2;
        new_dp = new_np + new_nn + 1;
        scratch = new_dp + qn + 1;

        if (cnt != 0)
        {
            cy = mpn_lshift(new_np, A + An - new_nn, new_nn, cnt);
            new_np[new_nn] = cy;
            new_nn += (cy != 0);
            mpn_lshift(new_dp, B + Bn - (qn + 1), qn + 1, cnt);
            new_dp[0] |= B[Bn - (qn + 1) - 1] >> (FLINT_BITS - cnt);
        }
        else
        {
            cy = 0;
            flint_mpn_copyi(new_np, A + An - new_nn, new_nn);
            new_dp = (mp_ptr) B + Bn - (qn + 1);
        }

        if (qn + 1 >= FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF)
        {
            /* {tp, new_nn - qn} with floor <= tp <= floor + 1 */
            divapprox(tp, new_np, new_nn, new_dp, qn + 1);
            qh = (cy != 0) ? tp[qn + 1] : tp[qn];
        }
        else
        {
            if (qn + 1 == 2)
            {
                qh = DIV_QR_2_NORM(tp, new_np, new_nn, new_dp);
            }
            else
            {
                dinv = flint_mpn_preinv1(new_dp[qn], new_dp[qn - 1]);
                if (qn + 1 < FLINT_MPN_DIVAPPR_DC_CUTOFF)
                    qh = divapprox_basecase_preinv1(tp, new_np, new_nn, new_dp, qn + 1, dinv);
                else
                    qh = divapprox_preinv1(tp, new_np, new_nn, new_dp, qn + 1, dinv, scratch);
            }
        }

        if (cy == 0)
        {
            tp[qn] = qh;
        }
        else if (FLINT_UNLIKELY(qh != 0))
        {
            /* the quotient is close to B^(qn + 1) and the approximation
               returned B^(qn + 1) */
            flint_mpn_store(tp, qn + 1, UWORD_MAX);
        }

        flint_mpn_copyi(Q, tp + 1, qn);

        if (exact && FLINT_UNLIKELY(tp[0] <= 4))
        {
            mp_size_t rn;

            rp = TMP_ALLOC((Bn + qn) * sizeof(mp_limb_t));
            FLINT_ASSERT(Bn > qn);   /* short quotient */
            flint_mpn_mul(rp, B, Bn, tp + 1, qn);
            rn = Bn + qn;
            rn -= (rp[rn - 1] == 0);

            if (rn > An || mpn_cmp(A, rp, An) < 0)
                mpn_sub_1(Q, Q, qn, 1);
        }

        TMP_END;
    }
}


/* the algorithm used by _flint_mpn_tdiv_qr (and flint_mpn_divapprox) */
enum { TDIV_SMALL, TDIV_NEWTON, TDIV_UNBALANCED, TDIV_PREINVN, TDIV_DIVCONQUER };

FLINT_FORCE_INLINE int
tdiv_method(mp_size_t An, mp_size_t Bn, int want_r)
{
    mp_size_t n = An - Bn + 1;

    if (FLINT_MPN_DIV_SMALL_SHAPE(An, Bn, want_r))
        return TDIV_SMALL;
    if ((Bn >= FLINT_MPN_TDIV_QR_NEWTON_CUTOFF && n >= FLINT_MPN_TDIV_QR_NEWTON_CUTOFF)
        || (Bn >= FLINT_MPN_TDIV_QR_NEWTON_LONG_CUTOFF && n >= 2 * Bn))
        return (An >= 3 * Bn) ? TDIV_UNBALANCED : TDIV_NEWTON;
    if ((Bn >= FLINT_MPN_TDIV_QR_UNBALANCED4_CUTOFF && An >= 4 * Bn)
          || (Bn >= FLINT_MPN_TDIV_QR_UNBALANCED3_CUTOFF && An >= 3 * Bn))
        return TDIV_UNBALANCED;
    if ((Bn >= 32 && An >= 4 * Bn) || (Bn >= 4 && An >= 32 * Bn))
        return TDIV_PREINVN;
    return TDIV_DIVCONQUER;
}

/*
    Approximate quotient floor(A / B) <= Q <= floor(A / B) + 1 with
    n = An - Bn + 1 limbs.
*/
/* approximate quotient by Karp-Markstein division with a guard limb */
static void
divapprox_newton(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn,
    mp_size_t f)
{
    /* as in _flint_mpn_tdiv_qr_newton, the integer quotient is {q, n}
       (with q[n] = 0) and the error is at most 4 units of the guard limb
       q[-1] (only 4/B at the integer scale). Adding 5 units of q[-1] and
       truncating gives floor(A/B) or floor(A/B) + 1. */
    /* with f fraction limbs: A B^f / B, the zero limbs not formed (the
       precision of _mp_real_div_newton is independent of An) */
    mp_size_t n = An + f - Bn + 1, n2 = n + 1;
    mp_ptr U, q;
    TMP_INIT;

    TMP_START;
    U = TMP_ALLOC((n2 + 2) * sizeof(mp_limb_t));
    _mp_real_div_newton(U, A, An, B, Bn, n2);
    q = U + n2 + 2 - (n + 1);
    mpn_add_1(q - 1, q - 1, n + 2, 5);
    if (FLINT_UNLIKELY(q[n] != 0))
        flint_mpn_store(Q, n, UWORD_MAX);   /* floor(A/B) = B^n - 1 */
    else
        flint_mpn_copyi(Q, q, n);
    TMP_END;
}

FLINT_FORCE_INLINE void
divapprox(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    mp_size_t n = An - Bn + 1;
    int method;

    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    /* the algorithm of the exact quotient, except that the Newton division
       and the divide and conquer division skip the correction steps, and
       that Karp-Markstein division with a guard limb is used from
       FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF instead of the divide and conquer
       division (unless the shape calls for the block or preinvn division),
       so that this is never slower than flint_mpn_tdiv_q */
    method = tdiv_method(An, Bn, 0);

    if (method == TDIV_SMALL)
    {
        tdiv_qr_small(Q, NULL, A, An, B, Bn);
    }
    else if (method == TDIV_UNBALANCED)
    {
        _flint_mpn_tdiv_qr_unbalanced(Q, NULL, A, An, B, Bn);
    }
    else if (method == TDIV_PREINVN)
    {
        _flint_mpn_tdiv_qr_preinvn(Q, NULL, A, An, B, Bn);
    }
    else if (method == TDIV_NEWTON ||
        (Bn >= FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF && n >= FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF))
    {
        divapprox_newton(Q, A, An, B, Bn, 0);
    }
    else
    {
        tdiv_q_divconquer(Q, A, An, B, Bn, 0);
    }
}

/*
    floor(A B^f / B) or floor(A B^f / B) + 1 with qn = An + f - Bn + 1
    limbs: f fraction limbs for f > 0 (without forming the zero limbs,
    except for long quotients outside the Newton range), only the top qn
    limbs of the quotient for f < 0.
*/
void
flint_mpn_divapprox_fraction(mp_ptr Q, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn, mp_size_t f)
{
    mp_size_t qn = An + f - Bn + 1;
    int method;

    FLINT_ASSERT(qn >= 1);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    if (f <= 0)
    {
        /* floor(A / (b B^-f)) = floor(floor(A / B^-f) / b) */
        divapprox(Q, A - f, An + f, B, Bn);
        return;
    }

    /* the algorithm of flint_mpn_divapprox on the padded numerator: the
       register-based, blockwise (unbalanced / preinvn) and divide and
       conquer divisions on the padded numerator, the short and Newton
       divisions without it */
    method = tdiv_method(An + f, Bn, 0);

    if (method != TDIV_SMALL && method != TDIV_UNBALANCED && method != TDIV_PREINVN &&
        (method == TDIV_NEWTON ||
        (Bn >= FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF && qn >= FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF)))
    {
        divapprox_newton(Q, A, An, B, Bn, f);
    }
    else if (method == TDIV_DIVCONQUER && Bn >= 2 &&
        qn + 1 >= FLINT_MPN_DIVAPPROX_SHORT_CUTOFF && qn <= Bn + 2)
    {
        tdiv_q_short(Q, A, An, B, Bn, f, 0);
    }
    else
    {
        /* short or long quotients: form A B^f */
        mp_ptr T;
        TMP_INIT;

        TMP_START;
        T = TMP_ALLOC((An + f) * sizeof(mp_limb_t));
        flint_mpn_zero(T, f);
        flint_mpn_copyi(T + f, A, An);
        divapprox(Q, T, An + f, B, Bn);
        TMP_END;
    }
}

/* R may be NULL */
void
_flint_mpn_tdiv_qr(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An,
    mp_srcptr B, mp_size_t Bn)
{
    FLINT_ASSERT(An >= Bn);
    FLINT_ASSERT(Bn >= 1);
    FLINT_ASSERT(B[Bn - 1] != 0);

    switch (tdiv_method(An, Bn, R != NULL))
    {
        case TDIV_SMALL:
            tdiv_qr_small(Q, R, A, An, B, Bn);
            break;
        case TDIV_NEWTON:
            _flint_mpn_tdiv_qr_newton(Q, R, A, An, B, Bn);
            break;
        case TDIV_UNBALANCED:
            _flint_mpn_tdiv_qr_unbalanced(Q, R, A, An, B, Bn);
            break;
        case TDIV_PREINVN:
            _flint_mpn_tdiv_qr_preinvn(Q, R, A, An, B, Bn);
            break;
        default:
            if (R != NULL)
                tdiv_qr_divconquer(Q, R, A, An, B, Bn);
            else
                tdiv_q_divconquer(Q, A, An, B, Bn, 1);
    }
}

/* exported versions *********************************************************/

void
_flint_mpn_tdiv_qr_small(mp_ptr qp, mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    tdiv_qr_small(qp, rp, ap, an, bp, bn);
}

mp_limb_t
_flint_mpn_divrem_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    return divrem_basecase_preinv1(qp, np, nn, dp, dn, dinv);
}

mp_limb_t
_flint_mpn_divapprox_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    return divapprox_basecase_preinv1(qp, np, nn, dp, dn, dinv);
}

mp_limb_t
_flint_mpn_div_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
    return div_basecase_preinv1(qp, np, nn, dp, dn, dinv);
}

mp_limb_t
_flint_mpn_divrem_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
{
    return divrem_n_divconquer_preinv1(qp, np, dp, n, dinv, tp);
}

mp_limb_t
_flint_mpn_divapprox_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
{
    return divapprox_n_divconquer_preinv1(qp, np, dp, n, dinv, tp);
}

mp_limb_t
_flint_mpn_divrem_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
{
    return divrem_preinv1(qp, np, nn, dp, dn, dinv, tp);
}

mp_limb_t
_flint_mpn_divapprox_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
{
    return divapprox_preinv1(qp, np, nn, dp, dn, dinv, tp);
}

void
_flint_mpn_tdiv_qr_divconquer(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    tdiv_qr_divconquer(Q, R, A, An, B, Bn);
}

void
_flint_mpn_tdiv_q_divconquer(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    tdiv_q_divconquer(Q, A, An, B, Bn, 1);
}

void
flint_mpn_divapprox(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)
{
    divapprox(Q, A, An, B, Bn);
}

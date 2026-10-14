/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Table-driven kernels for one to three fraction limbs (64-bit only).

    These evaluate the elementary functions on a fixed-point argument v
    in [0, 1) at n = 1, 2 or 3 fraction limbs, entirely in registers: the
    top bits of v index precomputed tables (small_tables.c), which leave
    a remainder below 2^-14 (n = 1, 2) resp. 2^-21 (n = 3) for a short
    polynomial, and the pieces are combined with truncated
    multiplications.  The one-limb kernels read the top limbs of the same
    table entries (floors at three limbs, hence floors at one).  Their contract
    is that of the mp_real kernels (impl.h), with smaller error bounds:
    the results are truncations within MP_REAL_SMALL_*_ERR ulps of B^-n.

    The cores are FLINT_FORCE_INLINE on scalar limbs so that callers
    with a fixed size (the nfloat fast paths) keep everything in
    registers; the mp_real kernels reach them through the array wrappers
    at the end.

    Notation: a k-limb fraction (a_{k-1}, ..., a_0) is
    sum a_i B^(i - k); "mulhi" keeps the top limbs of a product of
    fractions, truncated (an error of a few ulps of the output).
*/

#ifndef MP_REAL_SMALL_H
#define MP_REAL_SMALL_H

#include "longlong.h"
#include "small_tables.h"
#include "hand_mulhi.inc"

#if FLINT_BITS == 64

/* multiplication helpers *****************************************************/

/* All products below truncate (the result is at most the exact product).
   In ulps of the output: _mp_real_mulhi_2x1 < 1, _mp_real_mulhi_2x2_sloppy
   and _mp_real_sqrhi_2x2 < 3, _mp_real_mulhi_3x2 < 3, the 3x3 products
   below < 5, the 3x1 product < 1. The error bounds of the kernels
   (MP_REAL_SMALL_*_ERR) add these, the table truncations (every entry a
   floor, below one ulp at either size) and the propagated errors
   (factors below 3), with some slack; the dropped series terms are below
   1/2 ulp. */

/* top 3 limbs of the product of two 3-limb fractions: the partial
   products of weight B^-1 .. B^-3 (the high words of the B^-4 column),
   below 5 ulps (three low words and the two products of the B^-5
   column dropped) */
FLINT_FORCE_INLINE void
_mp_real_small_mulhi_3x3(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0, ulong b2, ulong b1, ulong b0)
{
    ulong t2, t1, x2, x1, x0, h1, h0, ph, pl, qh, ql;

    /* the B^-4 column (high words) and the B^-3 column summed apart from
       the top product, shortening the carry chains */
    h0 = n_mulhi(a2, b0);
    add_ssaaaa(h1, h0, UWORD(0), h0, UWORD(0), n_mulhi(a1, b1));
    add_ssaaaa(h1, h0, h1, h0, UWORD(0), n_mulhi(a0, b2));
    umul_ppmm(ph, pl, a2, b1);
    umul_ppmm(qh, ql, a1, b2);
    add_sssaaaaaa(x2, x1, x0, UWORD(0), ph, pl, UWORD(0), qh, ql);
    add_sssaaaaaa(x2, x1, x0, x2, x1, x0, UWORD(0), h1, h0);
    umul_ppmm(t2, t1, a2, b2);
    add_ssaaaa(*r2, *r1, t2, t1, x2, x1);
    *r0 = x0;
}

/* top 3 limbs of the square of a 3-limb fraction, below 5 ulps */
FLINT_FORCE_INLINE void
_mp_real_small_sqrhi_3x3(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0)
{
    ulong t2, t1, x2, x1, x0, h1, h0, ph, pl, c;

    c = n_mulhi(a2, a0);
    h1 = c >> (FLINT_BITS - 1);
    h0 = c << 1;
    add_ssaaaa(h1, h0, h1, h0, UWORD(0), n_mulhi(a1, a1));
    umul_ppmm(ph, pl, a2, a1);
    x2 = ph >> (FLINT_BITS - 1);
    x1 = (ph << 1) | (pl >> (FLINT_BITS - 1));
    x0 = pl << 1;
    add_sssaaaaaa(x2, x1, x0, x2, x1, x0, UWORD(0), h1, h0);
    umul_ppmm(t2, t1, a2, a2);
    add_ssaaaa(*r2, *r1, t2, t1, x2, x1);
    *r0 = x0;
}

/* top 3 limbs of a 3-limb fraction times a one-limb fraction b, below
   1 ulp */
FLINT_FORCE_INLINE void
_mp_real_small_mulhi_3x1(ulong * r2, ulong * r1, ulong * r0,
    ulong a2, ulong a1, ulong a0, ulong b)
{
    ulong t2, t1, ph, pl;

    umul_ppmm(t2, t1, a2, b);
    umul_ppmm(ph, pl, a1, b);
    add_sssaaaaaa(t2, t1, pl, t2, t1, UWORD(0), UWORD(0), ph, pl);
    ph = n_mulhi(a0, b);
    add_sssaaaaaa(t2, t1, pl, t2, t1, pl, UWORD(0), UWORD(0), ph);

    *r2 = t2;
    *r1 = t1;
    *r0 = pl;
}

/* exp ***********************************************************************/

/* 1/k! as fractions */
#define MP_REAL_SMALL_EXP_C3_2 UWORD(0x2aaaaaaaaaaaaaaa)   /* 1/6: B^-1 .. B^-3 */
#define MP_REAL_SMALL_EXP_C3_1 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_EXP_C3_0 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_EXP_C4_1 UWORD(0x0aaaaaaaaaaaaaaa)   /* 1/24 */
#define MP_REAL_SMALL_EXP_C4_0 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_EXP_C5_1 UWORD(0x0222222222222222)   /* 1/120 */
#define MP_REAL_SMALL_EXP_C5_0 UWORD(0x2222222222222222)
#define MP_REAL_SMALL_EXP_C6_1 UWORD(0x005b05b05b05b05b)   /* 1/720 */
#define MP_REAL_SMALL_EXP_C6_0 UWORD(0x05b05b05b05b05b0)
#define MP_REAL_SMALL_EXP_C7 UWORD(0x000d00d00d00d00d)     /* 1/5040 */
#define MP_REAL_SMALL_EXP_C8 UWORD(0x0001a01a01a01a01)     /* 1/40320 */

/* error bounds in ulps of B^-n (the results are below the exact values):
   at two limbs f within 4.5 ulps (r^2 3 ulps, times P < 1/2, and the
   product), c = b + f + b f within 8.6 (E2 one ulp, the product 3), y
   within A 8.6 + 3 + 1 < 27.1 with A = E1[j1] < e^(127/128) < 2.7; at
   three limbs f within 7.5, c within 13.5 and 19.6 after the two table
   levels, y within 2.7 * 19.6 + 5 + 1.01 < 59 */
#define MP_REAL_SMALL_EXP_ERR_2 32
#define MP_REAL_SMALL_EXP_ERR_3 64

/* one fraction limb: the polynomial within 1.7 ulps, c = b + f + b f
   within 3.8, y = A (1 + c) within 1.01 + 2.7 * 3.8 + 1 < 12.3 */
#define MP_REAL_SMALL_EXP_ERR_1 16

/*
    (y1, y0) = exp(v) at one fraction limb (y1 the integral limb, 1 or 2)
    for v in [0, 1): as at two limbs,

        v = j1/2^7 + j2/2^14 + r,  r < 2^-14,
        exp(v) = E1[j1] (1 + E2[j2]) (1 + f),  f = r + r^2 (1/2 + r/6 + r^2/24)

    (the next term, r^5/120, below 2^-76), every product one high word.
*/
FLINT_FORCE_INLINE void
_mp_real_small_exp_1(ulong * y1, ulong * y0, ulong v)
{
    ulong j1, j2, R, q, p, f, b, c, a1, a0, u1, u0;
    const ulong * A;

    j1 = v >> 57;
    j2 = (v >> 50) & 127;
    R = v & ((UWORD(1) << 50) - 1);

    q = n_mulhi(R, R);
    p = MP_REAL_SMALL_EXP_C3_2 + n_mulhi(R, MP_REAL_SMALL_EXP_C4_1);
    p = (UWORD(1) << 63) + n_mulhi(R, p);
    f = R + n_mulhi(q, p);

    /* c = (1 + b)(1 + f) - 1, b = E2[j2] */
    b = _mp_real_small_exp_e2[3 * j2 + 2];
    c = b + f + n_mulhi(b, f);

    /* y = A (1 + c) = a1 + a0 + a1 c + a0 c, a1 in {1, 2} */
    A = _mp_real_small_exp_e1 + 4 * j1;
    a1 = A[3];
    a0 = A[2];
    u1 = (a1 == 2) ? (c >> 63) : 0;
    u0 = (a1 == 2) ? (c << 1) : c;
    add_ssaaaa(a1, a0, a1, a0, u1, u0);
    add_ssaaaa(*y1, *y0, a1, a0, UWORD(0), n_mulhi(A[2], c));
}

/*
    (y2, y1, y0) = exp(v) at two fraction limbs (y2 the integral limb,
    1 or 2) for v = (v1, v0) in [0, 1):

        v = j1/2^7 + j2/2^14 + r,  r < 2^-14,
        exp(v) = E1[j1] (1 + E2[j2]) (1 + f),  f = exp(r) - 1,

    f = r + r^2 P(r) with P = 1/2 + r/6 + ... + r^6/8! (the next term,
    r^9/9!, below 2^-144), Horner from the top at one limb (on the top
    64 bits of r 2^14) while the terms allow it, then at two.
*/
FLINT_FORCE_INLINE void
_mp_real_small_exp_2(ulong * y2, ulong * y1, ulong * y0, ulong v1, ulong v0)
{
    ulong j1, j2, r1, r0, rho, g, h1, h0, q1, q0, p1, p0, s1, s0, f1, f0;
    ulong b1, b0, c1, c0, t1, t0, a2, a1, a0, u2, u1, u0;
    const ulong * A;
    const ulong * E;

    j1 = v1 >> 57;
    j2 = (v1 >> 50) & 127;
    r1 = v1 & ((UWORD(1) << 50) - 1);
    r0 = v0;
    rho = (r1 << 14) | (r0 >> 50);          /* floor(r 2^78) */

    /* G = 1/120 + r/720 + r^2/5040 + r^3/8!, at one limb (B^-1) */
    g = MP_REAL_SMALL_EXP_C7 + (n_mulhi(rho, MP_REAL_SMALL_EXP_C8) >> 14);
    g = MP_REAL_SMALL_EXP_C6_1 + (n_mulhi(rho, g) >> 14);
    g = MP_REAL_SMALL_EXP_C5_1 + (n_mulhi(rho, g) >> 14);

    /* H = 1/24 + r G, Q = 1/6 + r H, P = 1/2 + r Q at two limbs */
    _mp_real_mulhi_2x1(&h1, &h0, r1, r0, g);
    add_ssaaaa(h1, h0, h1, h0, MP_REAL_SMALL_EXP_C4_1, MP_REAL_SMALL_EXP_C4_0);
    _mp_real_mulhi_2x2_sloppy(&q1, &q0, r1, r0, h1, h0);
    add_ssaaaa(q1, q0, q1, q0, MP_REAL_SMALL_EXP_C3_2, MP_REAL_SMALL_EXP_C3_1);
    _mp_real_mulhi_2x2_sloppy(&p1, &p0, r1, r0, q1, q0);
    p1 += UWORD(1) << 63;

    /* f = r + r^2 P */
    _mp_real_sqrhi_2x2(&s1, &s0, r1, r0);
    _mp_real_mulhi_2x2_sloppy(&s1, &s0, s1, s0, p1, p0);
    add_ssaaaa(f1, f0, r1, r0, s1, s0);

    /* c = (1 + b)(1 + f) - 1 = b + f + b f, b = E2[j2] */
    E = _mp_real_small_exp_e2 + 3 * j2;
    b1 = E[2];
    b0 = E[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, b1, b0, f1, f0);
    add_ssaaaa(c1, c0, b1, b0, f1, f0);
    add_ssaaaa(c1, c0, c1, c0, t1, t0);

    /* y = A (1 + c) = A + a2 c + (a1, a0) c, A = E1[j1] = a2 + (a1, a0),
       a2 in {1, 2} */
    A = _mp_real_small_exp_e1 + 4 * j1;
    a2 = A[3];
    a1 = A[2];
    a0 = A[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, a1, a0, c1, c0);
    u2 = 0;
    u1 = c1;
    u0 = c0;
    if (a2 == 2)
    {
        u2 = c1 >> 63;
        u1 = (c1 << 1) | (c0 >> 63);
        u0 = c0 << 1;
    }
    add_sssaaaaaa(a2, a1, a0, a2, a1, a0, u2, u1, u0);
    add_sssaaaaaa(*y2, *y1, *y0, a2, a1, a0, UWORD(0), t1, t0);
}

/*
    (y3, y2, y1, y0) = exp(v) at three fraction limbs for
    v = (v2, v1, v0) in [0, 1):

        v = j1/2^7 + j2/2^14 + j3/2^21 + r,  r < 2^-21,
        exp(v) = E1[j1] (1 + E2[j2]) (1 + E3[j3]) (1 + f),

    f = r + r^2 P(r), P = 1/2 + r/6 + ... + r^6/8! (the next term below
    2^-207): one limb for the innermost step, two limbs while the
    remaining powers of r allow it, then three.
*/
FLINT_FORCE_INLINE void
_mp_real_small_exp_3(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong v2, ulong v1, ulong v0)
{
    ulong j1, j2, j3, r2, r1, r0, rho, e, x1, x0, q2, q1, q0, s2, s1, s0;
    ulong f2, f1, f0, b2, b1, b0, c2, c1, c0, t2, t1, t0, a3, a2, a1, a0;
    ulong u3, u2, u1, u0;
    const ulong * A;
    const ulong * E;

    j1 = v2 >> 57;
    j2 = (v2 >> 50) & 127;
    j3 = (v2 >> 43) & 127;
    r2 = v2 & ((UWORD(1) << 43) - 1);
    r1 = v1;
    r0 = v0;
    rho = (r2 << 21) | (r1 >> 43);          /* floor(r 2^85) */

    /* 1/5040 + r/8! at one limb, then 1/720 + r (..), 1/120 + r (..),
       1/24 + r (..) at two limbs (on the top two limbs of r) */
    e = MP_REAL_SMALL_EXP_C7 + (n_mulhi(rho, MP_REAL_SMALL_EXP_C8) >> 21);
    _mp_real_mulhi_2x1(&x1, &x0, r2, r1, e);
    add_ssaaaa(x1, x0, x1, x0, MP_REAL_SMALL_EXP_C6_1, MP_REAL_SMALL_EXP_C6_0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, r2, r1, x1, x0);
    add_ssaaaa(x1, x0, x1, x0, MP_REAL_SMALL_EXP_C5_1, MP_REAL_SMALL_EXP_C5_0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, r2, r1, x1, x0);
    add_ssaaaa(x1, x0, x1, x0, MP_REAL_SMALL_EXP_C4_1, MP_REAL_SMALL_EXP_C4_0);

    /* Q = 1/6 + r H, P = 1/2 + r Q at three limbs */
    _mp_real_mulhi_3x2(&q2, &q1, &q0, r2, r1, r0, x1, x0);
    add_sssaaaaaa(q2, q1, q0, q2, q1, q0, MP_REAL_SMALL_EXP_C3_2, MP_REAL_SMALL_EXP_C3_1, MP_REAL_SMALL_EXP_C3_0);
    _mp_real_small_mulhi_3x3(&q2, &q1, &q0, r2, r1, r0, q2, q1, q0);
    q2 += UWORD(1) << 63;

    /* f = r + r^2 P */
    _mp_real_small_sqrhi_3x3(&s2, &s1, &s0, r2, r1, r0);
    _mp_real_small_mulhi_3x3(&s2, &s1, &s0, s2, s1, s0, q2, q1, q0);
    add_sssaaaaaa(f2, f1, f0, r2, r1, r0, s2, s1, s0);

    /* c = (1 + E3)(1 + f) - 1, then (1 + E2)(1 + c) - 1 */
    E = _mp_real_small_exp_e3 + 3 * j3;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, f2, f1, f0);
    add_sssaaaaaa(c2, c1, c0, b2, b1, b0, f2, f1, f0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, t2, t1, t0);

    E = _mp_real_small_exp_e2 + 3 * j2;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, c2, c1, c0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, b2, b1, b0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, t2, t1, t0);

    /* y = A (1 + c), A = E1[j1] */
    A = _mp_real_small_exp_e1 + 4 * j1;
    a3 = A[3]; a2 = A[2]; a1 = A[1]; a0 = A[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, a2, a1, a0, c2, c1, c0);
    u3 = 0; u2 = c2; u1 = c1; u0 = c0;
    if (a3 == 2)
    {
        u3 = c2 >> 63;
        u2 = (c2 << 1) | (c1 >> 63);
        u1 = (c1 << 1) | (c0 >> 63);
        u0 = c0 << 1;
    }
    add_ssssaaaaaaaa(a3, a2, a1, a0, a3, a2, a1, a0, u3, u2, u1, u0);
    add_ssssaaaaaaaa(*y3, *y2, *y1, *y0, a3, a2, a1, a0, UWORD(0), t2, t1, t0);
}

/* array wrappers: (y, n + 1) = exp(v), v at n = 1, 2 resp. 3 fraction
   limbs */
FLINT_FORCE_INLINE void
_mp_real_small_exp_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_exp_1(y + 1, y, v[0]);
    *err = MP_REAL_SMALL_EXP_ERR_1;
}

FLINT_FORCE_INLINE void
_mp_real_small_exp_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_exp_2(y + 2, y + 1, y, v[1], v[0]);
    *err = MP_REAL_SMALL_EXP_ERR_2;
}

FLINT_FORCE_INLINE void
_mp_real_small_exp_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_exp_3(y + 3, y + 2, y + 1, y, v[2], v[1], v[0]);
    *err = MP_REAL_SMALL_EXP_ERR_3;
}

/*
    exp(t) and exp(-t) together, for the hyperbolic functions: with the
    same splitting of t, the even and odd parts of the polynomial,

        exp(+-r) - 1 = +-r + s (Pe +- r Po),  s = r^2,
        Pe = 1/2 + s/24 + s^2/720 + s^3/8!,  Po = 1/6 + s/120 + s^2/5040,

    serve both signs (the dropped terms as for exp alone), and
    exp(-t) = N1[j1] (1 - N2[j2]) (1 - g) with g = 1 - exp(-r) (and N3 at
    three limbs) from the tables of exp(-j/2^7) and 1 - exp(-j/2^14),
    1 - exp(-j/2^21). Both outputs are within the bounds of the kernels
    for exp: exp(t) as there (the polynomial at most 7 ulps off where
    it enters multiplied by s), exp(-t) <= 1 within 12.6 resp. 25.6 ulps
    (g within 4.5 resp. 7.5, then one floor and one product per table
    level, the errors not amplified by factors above 1). The shared
    work makes the pair about 1.5 times the cost of one exp kernel,
    against 2 for two kernels; the two halves are independent after
    the polynomial, for instruction-level parallelism.
*/

/* (y2, y1, y0) = exp(t), (z2, z1, z0) = exp(-t) (the integral limbs first)
   for t = (v1, v0) in [0, 1) */
FLINT_FORCE_INLINE void
_mp_real_small_exp_pm_2(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0,
    ulong v1, ulong v0)
{
    ulong j1, j2, r1, r0, s1, s0, sigma, ge, go, h1, h0, o1, o0, e1, e0, p1, p0, m1, m0, w1, w0;
    ulong fp1, fp0, gm1, gm0, b1, b0, c1, c0, t1, t0, a2, a1, a0, u2, u1, u0;
    const ulong * A;
    const ulong * E;

    j1 = v1 >> 57;
    j2 = (v1 >> 50) & 127;
    r1 = v1 & ((UWORD(1) << 50) - 1);
    r0 = v0;

    _mp_real_sqrhi_2x2(&s1, &s0, r1, r0);          /* s = r^2 < 2^-28 */
    sigma = (s1 << 28) | (s0 >> 36);               /* floor(s 2^92) */

    /* 1/720 + s/8!, 1/120 + s/5040 at one limb; then 1/24 + s (..) and
       Po = 1/6 + s (..) at two, Pe = 1/2 + s (..) */
    ge = MP_REAL_SMALL_EXP_C6_1 + (n_mulhi(sigma, MP_REAL_SMALL_EXP_C8) >> 28);
    go = MP_REAL_SMALL_EXP_C5_1 + (n_mulhi(sigma, MP_REAL_SMALL_EXP_C7) >> 28);
    _mp_real_mulhi_2x1(&h1, &h0, s1, s0, ge);
    add_ssaaaa(h1, h0, h1, h0, MP_REAL_SMALL_EXP_C4_1, MP_REAL_SMALL_EXP_C4_0);
    _mp_real_mulhi_2x1(&o1, &o0, s1, s0, go);
    add_ssaaaa(o1, o0, o1, o0, MP_REAL_SMALL_EXP_C3_2, MP_REAL_SMALL_EXP_C3_1);
    _mp_real_mulhi_2x2_sloppy(&e1, &e0, s1, s0, h1, h0);
    e1 += UWORD(1) << (FLINT_BITS - 1);
    _mp_real_mulhi_2x2_sloppy(&o1, &o0, r1, r0, o1, o0);

    /* exp(r) - 1 = r + s (Pe + r Po), 1 - exp(-r) = r - s (Pe - r Po) */
    add_ssaaaa(p1, p0, e1, e0, o1, o0);
    sub_ddmmss(m1, m0, e1, e0, o1, o0);
    _mp_real_mulhi_2x2_sloppy(&w1, &w0, s1, s0, p1, p0);
    add_ssaaaa(fp1, fp0, r1, r0, w1, w0);
    _mp_real_mulhi_2x2_sloppy(&w1, &w0, s1, s0, m1, m0);
    sub_ddmmss(gm1, gm0, r1, r0, w1, w0);

    /* exp(t) = E1 (1 + c), c = b + f + b f, b = E2 */
    E = _mp_real_small_exp_e2 + 3 * j2;
    b1 = E[2];
    b0 = E[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, b1, b0, fp1, fp0);
    add_ssaaaa(c1, c0, b1, b0, fp1, fp0);
    add_ssaaaa(c1, c0, c1, c0, t1, t0);
    A = _mp_real_small_exp_e1 + 4 * j1;
    a2 = A[3];
    a1 = A[2];
    a0 = A[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, a1, a0, c1, c0);
    u2 = 0;
    u1 = c1;
    u0 = c0;
    if (a2 == 2)
    {
        u2 = c1 >> (FLINT_BITS - 1);
        u1 = (c1 << 1) | (c0 >> (FLINT_BITS - 1));
        u0 = c0 << 1;
    }
    add_sssaaaaaa(a2, a1, a0, a2, a1, a0, u2, u1, u0);
    add_sssaaaaaa(*y2, *y1, *y0, a2, a1, a0, UWORD(0), t1, t0);

    /* exp(-t) = N1 (1 - c'), c' = b + g - b g, b = N2 */
    E = _mp_real_small_exp_n2 + 3 * j2;
    b1 = E[2];
    b0 = E[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, b1, b0, gm1, gm0);
    add_ssaaaa(c1, c0, b1, b0, gm1, gm0);
    sub_ddmmss(c1, c0, c1, c0, t1, t0);
    A = _mp_real_small_exp_n1 + 4 * j1;
    a2 = A[3];
    a1 = A[2];
    a0 = A[1];
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, a1, a0, c1, c0);
    if (a2 != 0)    /* N1[0] = 1 */
        add_ssaaaa(t1, t0, t1, t0, c1, c0);
    sub_dddmmmsss(*z2, *z1, *z0, a2, a1, a0, UWORD(0), t1, t0);
}

/* the same at three limbs: exp(t) = (y3, y2, y1, y0), exp(-t) = (z3, z2, z1, z0) */
FLINT_FORCE_INLINE void
_mp_real_small_exp_pm_3(ulong * y3, ulong * y2, ulong * y1, ulong * y0,
    ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong v2, ulong v1, ulong v0)
{
    ulong j1, j2, j3, r2, r1, r0, s2, s1, s0, sigma, m, ge1, ge0, go1, go0, h1, h0;
    ulong o2, o1, o0, e2, e1, e0, p2, p1, p0, n2, n1, n0, w2, w1, w0;
    ulong f2, f1, f0, g2, g1, g0, b2, b1, b0, c2, c1, c0, t2, t1, t0;
    ulong a3, a2, a1, a0, u3, u2, u1, u0;
    const ulong * A;
    const ulong * E;

    j1 = v2 >> 57;
    j2 = (v2 >> 50) & 127;
    j3 = (v2 >> 43) & 127;
    r2 = v2 & ((UWORD(1) << 43) - 1);
    r1 = v1;
    r0 = v0;

    _mp_real_small_sqrhi_3x3(&s2, &s1, &s0, r2, r1, r0);  /* s = r^2 < 2^-42 */
    sigma = (s2 << 42) | (s1 >> 22);                      /* floor(s 2^106) */

    /* 1/720 + s/8!, 1/120 + s/5040 at two limbs (the s terms from one
       limb at the scale 2^-106: their errors enter multiplied by s^3
       resp. r s^2) */
    m = n_mulhi(sigma, MP_REAL_SMALL_EXP_C8);
    add_ssaaaa(ge1, ge0, MP_REAL_SMALL_EXP_C6_1, MP_REAL_SMALL_EXP_C6_0, m >> 42, m << 22);
    m = n_mulhi(sigma, MP_REAL_SMALL_EXP_C7);
    add_ssaaaa(go1, go0, MP_REAL_SMALL_EXP_C5_1, MP_REAL_SMALL_EXP_C5_0, m >> 42, m << 22);

    /* 1/24 + s (..) at two limbs, Pe = 1/2 + s (..) and Po = 1/6 + s (..)
       at three */
    _mp_real_mulhi_2x2_sloppy(&h1, &h0, s2, s1, ge1, ge0);
    add_ssaaaa(h1, h0, h1, h0, MP_REAL_SMALL_EXP_C4_1, MP_REAL_SMALL_EXP_C4_0);
    _mp_real_mulhi_3x2(&e2, &e1, &e0, s2, s1, s0, h1, h0);
    e2 += UWORD(1) << (FLINT_BITS - 1);
    _mp_real_mulhi_3x2(&o2, &o1, &o0, s2, s1, s0, go1, go0);
    add_sssaaaaaa(o2, o1, o0, o2, o1, o0, MP_REAL_SMALL_EXP_C3_2, MP_REAL_SMALL_EXP_C3_1, MP_REAL_SMALL_EXP_C3_0);
    _mp_real_small_mulhi_3x3(&o2, &o1, &o0, r2, r1, r0, o2, o1, o0);

    add_sssaaaaaa(p2, p1, p0, e2, e1, e0, o2, o1, o0);
    sub_dddmmmsss(n2, n1, n0, e2, e1, e0, o2, o1, o0);
    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, s2, s1, s0, p2, p1, p0);
    add_sssaaaaaa(f2, f1, f0, r2, r1, r0, w2, w1, w0);
    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, s2, s1, s0, n2, n1, n0);
    sub_dddmmmsss(g2, g1, g0, r2, r1, r0, w2, w1, w0);

    /* exp(t): c = (1 + E3)(1 + f) - 1, then (1 + E2)(1 + c) - 1, E1 (1 + c) */
    E = _mp_real_small_exp_e3 + 3 * j3;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, f2, f1, f0);
    add_sssaaaaaa(c2, c1, c0, b2, b1, b0, f2, f1, f0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, t2, t1, t0);
    E = _mp_real_small_exp_e2 + 3 * j2;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, c2, c1, c0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, b2, b1, b0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, t2, t1, t0);
    A = _mp_real_small_exp_e1 + 4 * j1;
    a3 = A[3]; a2 = A[2]; a1 = A[1]; a0 = A[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, a2, a1, a0, c2, c1, c0);
    u3 = 0; u2 = c2; u1 = c1; u0 = c0;
    if (a3 == 2)
    {
        u3 = c2 >> (FLINT_BITS - 1);
        u2 = (c2 << 1) | (c1 >> (FLINT_BITS - 1));
        u1 = (c1 << 1) | (c0 >> (FLINT_BITS - 1));
        u0 = c0 << 1;
    }
    add_ssssaaaaaaaa(a3, a2, a1, a0, a3, a2, a1, a0, u3, u2, u1, u0);
    add_ssssaaaaaaaa(*y3, *y2, *y1, *y0, a3, a2, a1, a0, UWORD(0), t2, t1, t0);

    /* exp(-t): c' = N3 + g - N3 g, then N2 + c' - N2 c', N1 (1 - c') */
    E = _mp_real_small_exp_n3 + 3 * j3;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, g2, g1, g0);
    add_sssaaaaaa(c2, c1, c0, b2, b1, b0, g2, g1, g0);
    sub_dddmmmsss(c2, c1, c0, c2, c1, c0, t2, t1, t0);
    E = _mp_real_small_exp_n2 + 3 * j2;
    b2 = E[2]; b1 = E[1]; b0 = E[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, b2, b1, b0, c2, c1, c0);
    add_sssaaaaaa(c2, c1, c0, c2, c1, c0, b2, b1, b0);
    sub_dddmmmsss(c2, c1, c0, c2, c1, c0, t2, t1, t0);
    A = _mp_real_small_exp_n1 + 4 * j1;
    a3 = A[3]; a2 = A[2]; a1 = A[1]; a0 = A[0];
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, a2, a1, a0, c2, c1, c0);
    if (a3 != 0)
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, c2, c1, c0);
    sub_ddddmmmmssss(*z3, *z2, *z1, *z0, a3, a2, a1, a0, UWORD(0), t2, t1, t0);
}

/* log1p *********************************************************************/

#define MP_REAL_SMALL_LOG_C2 UWORD(0x8000000000000000)     /* 1/2 */
#define MP_REAL_SMALL_LOG_C3 UWORD(0x5555555555555555)     /* 1/3 (every limb) */
#define MP_REAL_SMALL_LOG_C4 UWORD(0x4000000000000000)     /* 1/4 */
#define MP_REAL_SMALL_LOG_C5 UWORD(0x3333333333333333)     /* 1/5 (every limb) */
#define MP_REAL_SMALL_LOG_C6 UWORD(0x2aaaaaaaaaaaaaaa)     /* 1/6, top limb */

/* error bounds in ulps of B^-n: each level of the reduction leaves u an
   upper bound (by below one ulp per level, log1p' <= 1) and the sum of
   the floors L below the exact one (one ulp per level); u^2 P is below
   its value by under 4.5 ulps at two limbs (7.5 at three): the error lies
   in (-3, 7.5) at two limbs, (-4, 11.5) at three */
#define MP_REAL_SMALL_LOG_ERR_2 16
#define MP_REAL_SMALL_LOG_ERR_3 24

/* one level of the multiplicative reduction at two limbs: with
   j = floor(u 2^(s + 7)) (sh = 64 - s - 7 the shift of the top limb),
   u <- (1 + u)(1 - D[j]) - 1 = u - D - u D (u D truncated, so that u
   stays an upper bound; it remains >= 0) and acc += L[j] */
#define MP_REAL_SMALL_LOG_LEVEL_2(dtab, ltab, sh) \
    do { \
        ulong __j = u1 >> (sh), __d = (dtab)[__j], __p1, __p0; \
        const ulong * __L = (ltab) + 3 * __j; \
        _mp_real_mulhi_2x1(&__p1, &__p0, u1, u0, __d); \
        sub_ddmmss(u1, u0, u1, u0, __d, UWORD(0)); \
        sub_ddmmss(u1, u0, u1, u0, __p1, __p0); \
        add_ssaaaa(a1, a0, a1, a0, __L[2], __L[1]); \
    } while (0)

#define MP_REAL_SMALL_LOG_LEVEL_3(dtab, ltab, sh) \
    do { \
        ulong __j = u2 >> (sh), __d = (dtab)[__j], __p2, __p1, __p0; \
        const ulong * __L = (ltab) + 3 * __j; \
        _mp_real_small_mulhi_3x1(&__p2, &__p1, &__p0, u2, u1, u0, __d); \
        sub_dddmmmsss(u2, u1, u0, u2, u1, u0, __d, UWORD(0), UWORD(0)); \
        sub_dddmmmsss(u2, u1, u0, u2, u1, u0, __p2, __p1, __p0); \
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, __L[2], __L[1], __L[0]); \
    } while (0)

/* one fraction limb: the three floors of L below by under 3 ulps, u an
   upper bound by under 3 (log1p' <= 1), u^2 P below by under 1.5 and the
   dropped -u^4/4 below 2^-86: within (-3, 4.5) */
#define MP_REAL_SMALL_LOG_ERR_1 8

/*
    log(1 + v) at one fraction limb for v in [0, 1): the three levels of
    the two-limb kernel (u < 2^-21 (1 + 2^-16)), then
    log(1 + u) = u - u^2 (1/2 - u/3) (the next term, u^4/4, below 2^-86).
*/
FLINT_FORCE_INLINE ulong
_mp_real_small_log1p_1(ulong v)
{
    ulong u = v, a = 0, j, d, p;

#define MP_REAL_SMALL_LOG_LEVEL_1(dtab, ltab, sh) \
    do { \
        j = u >> (sh); \
        d = (dtab)[j]; \
        p = n_mulhi(u, d); \
        u = u - d - p; \
        a += (ltab)[3 * j + 2]; \
    } while (0)

    MP_REAL_SMALL_LOG_LEVEL_1(_mp_real_small_log_d1, _mp_real_small_log_l1, 57);
    MP_REAL_SMALL_LOG_LEVEL_1(_mp_real_small_log_d2, _mp_real_small_log_l2, 50);
    MP_REAL_SMALL_LOG_LEVEL_1(_mp_real_small_log_d3, _mp_real_small_log_l3, 43);

#undef MP_REAL_SMALL_LOG_LEVEL_1

    p = MP_REAL_SMALL_LOG_C2 - n_mulhi(u, MP_REAL_SMALL_LOG_C3);
    return a + u - n_mulhi(n_mulhi(u, u), p);
}

/*
    (y1, y0) = log(1 + v) at two fraction limbs for v = (v1, v0) in
    [0, 1): three levels of the reduction (1 + v) prod (1 - D_k[j_k]) =
    1 + u with 0 <= u < 2^-21 (1 + 2^-16), so that

        log(1 + v) = sum L_k[j_k] + log(1 + u),

    and log(1 + u) = u - u^2 P(u), P = 1/2 - u/3 + u^2/4 - u^3/5 + u^4/6
    (the next term, u^7/7, below 2^-149), Horner at one limb for the
    inner terms (on the top bits of u 2^21), then at two.
*/
FLINT_FORCE_INLINE void
_mp_real_small_log1p_2(ulong * y1, ulong * y0, ulong v1, ulong v0)
{
    ulong u1 = v1, u0 = v0, a1 = 0, a0 = 0, rho, r, q1, q0, p1, p0, s1, s0;

    MP_REAL_SMALL_LOG_LEVEL_2(_mp_real_small_log_d1, _mp_real_small_log_l1, 57);
    MP_REAL_SMALL_LOG_LEVEL_2(_mp_real_small_log_d2, _mp_real_small_log_l2, 50);
    MP_REAL_SMALL_LOG_LEVEL_2(_mp_real_small_log_d3, _mp_real_small_log_l3, 43);

    /* R = 1/4 - u/5 + u^2/6 at one limb (u may exceed 2^-21 slightly:
       one bit of headroom) */
    rho = (u1 << 20) | (u0 >> 44);          /* floor(u 2^84) */
    r = MP_REAL_SMALL_LOG_C5 - (n_mulhi(rho, MP_REAL_SMALL_LOG_C6) >> 20);
    r = MP_REAL_SMALL_LOG_C4 - (n_mulhi(rho, r) >> 20);

    /* Q = 1/3 - u R, P = 1/2 - u Q */
    _mp_real_mulhi_2x1(&q1, &q0, u1, u0, r);
    sub_ddmmss(q1, q0, MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3, q1, q0);
    _mp_real_mulhi_2x2_sloppy(&p1, &p0, u1, u0, q1, q0);
    sub_ddmmss(p1, p0, MP_REAL_SMALL_LOG_C2, UWORD(0), p1, p0);

    /* y = acc + u - u^2 P */
    _mp_real_sqrhi_2x2(&s1, &s0, u1, u0);
    _mp_real_mulhi_2x2_sloppy(&s1, &s0, s1, s0, p1, p0);
    add_ssaaaa(a1, a0, a1, a0, u1, u0);
    sub_ddmmss(*y1, *y0, a1, a0, s1, s0);
}

/* (y2, y1, y0) = log(1 + v) at three fraction limbs: four levels of the
   reduction (u < 2^-28 (1 + 2^-16)), then the same polynomial (the next
   term below 2^-198) */
FLINT_FORCE_INLINE void
_mp_real_small_log1p_3(ulong * y2, ulong * y1, ulong * y0, ulong v2, ulong v1, ulong v0)
{
    ulong u2 = v2, u1 = v1, u0 = v0, a2 = 0, a1 = 0, a0 = 0, rho, r;
    ulong x1, x0, q1, q0, p2, p1, p0, s2, s1, s0;

    MP_REAL_SMALL_LOG_LEVEL_3(_mp_real_small_log_d1, _mp_real_small_log_l1, 57);
    MP_REAL_SMALL_LOG_LEVEL_3(_mp_real_small_log_d2, _mp_real_small_log_l2, 50);
    MP_REAL_SMALL_LOG_LEVEL_3(_mp_real_small_log_d3, _mp_real_small_log_l3, 43);
    MP_REAL_SMALL_LOG_LEVEL_3(_mp_real_small_log_d4, _mp_real_small_log_l4, 36);

    /* 1/5 - u/6 at one limb, R = 1/4 - u (..), Q = 1/3 - u R at two
       (on the top two limbs of u), P = 1/2 - u Q at three */
    rho = (u2 << 27) | (u1 >> 37);          /* floor(u 2^91) */
    r = MP_REAL_SMALL_LOG_C5 - (n_mulhi(rho, MP_REAL_SMALL_LOG_C6) >> 27);
    _mp_real_mulhi_2x1(&x1, &x0, u2, u1, r);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_LOG_C4, UWORD(0), x1, x0);
    _mp_real_mulhi_2x2_sloppy(&q1, &q0, u2, u1, x1, x0);
    sub_ddmmss(q1, q0, MP_REAL_SMALL_LOG_C3, MP_REAL_SMALL_LOG_C3, q1, q0);
    _mp_real_mulhi_3x2(&p2, &p1, &p0, u2, u1, u0, q1, q0);
    sub_dddmmmsss(p2, p1, p0, MP_REAL_SMALL_LOG_C2, UWORD(0), UWORD(0), p2, p1, p0);

    /* y = acc + u - u^2 P */
    _mp_real_small_sqrhi_3x3(&s2, &s1, &s0, u2, u1, u0);
    _mp_real_small_mulhi_3x3(&s2, &s1, &s0, s2, s1, s0, p2, p1, p0);
    add_sssaaaaaa(a2, a1, a0, a2, a1, a0, u2, u1, u0);
    sub_dddmmmsss(*y2, *y1, *y0, a2, a1, a0, s2, s1, s0);
}

/* array wrappers: (y, n) = log1p(v), v at n = 2 resp. 3 fraction limbs */
FLINT_FORCE_INLINE void
_mp_real_small_log1p_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    y[0] = _mp_real_small_log1p_1(v[0]);
    *err = MP_REAL_SMALL_LOG_ERR_1;
}

FLINT_FORCE_INLINE void
_mp_real_small_log1p_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_log1p_2(y + 1, y, v[1], v[0]);
    *err = MP_REAL_SMALL_LOG_ERR_2;
}

FLINT_FORCE_INLINE void
_mp_real_small_log1p_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_log1p_3(y + 2, y + 1, y, v[2], v[1], v[0]);
    *err = MP_REAL_SMALL_LOG_ERR_3;
}

/* sin and cos *****************************************************************/

/* 1/k! (top limbs first) */
#define MP_REAL_SMALL_F6_2 UWORD(0x2aaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F6_1 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F6_0 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F24_2 UWORD(0x0aaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F24_1 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F24_0 UWORD(0xaaaaaaaaaaaaaaaa)
#define MP_REAL_SMALL_F120_2 UWORD(0x0222222222222222)
#define MP_REAL_SMALL_F120_1 UWORD(0x2222222222222222)
#define MP_REAL_SMALL_F720_2 UWORD(0x005b05b05b05b05b)
#define MP_REAL_SMALL_F720_1 UWORD(0x05b05b05b05b05b0)
#define MP_REAL_SMALL_F5040_2 UWORD(0x000d00d00d00d00d)
#define MP_REAL_SMALL_F5040_1 UWORD(0x00d00d00d00d00d0)
#define MP_REAL_SMALL_F8_2 UWORD(0x0001a01a01a01a01)     /* 1/8! */
#define MP_REAL_SMALL_F8_1 UWORD(0xa01a01a01a01a01a)
#define MP_REAL_SMALL_F9_2 UWORD(0x00002e3bc74aad8e)     /* 1/9! */
#define MP_REAL_SMALL_F9_1 UWORD(0x671f5583911ca002)
#define MP_REAL_SMALL_F10_2 UWORD(0x0000049f93edde27)    /* 1/10! */
#define MP_REAL_SMALL_F11_2 UWORD(0x0000006b99159fd5)    /* 1/11! */

/* error bounds in ulps of B^-n, for both outputs: the table pairs
   combined within 9.4 ulps at two limbs (13.3 at three: the floors, two
   products per output, and the errors of one table times the values of
   the other, |sin|, 1 - cos < 0.85), the polynomial parts within 3.6 and
   5.2 (5.9 and 8.3), and the final angle addition within 26 (40.3) */
#define MP_REAL_SMALL_SIN_COS_ERR_2 32
#define MP_REAL_SMALL_SIN_COS_ERR_3 48

/* angle addition on (sin, 1 - cos) pairs at two limbs:
   sin(a + b) = sa + sb - sa gb - ga sb,
   1 - cos(a + b) = ga + gb - ga gb + sa sb */
FLINT_FORCE_INLINE void
_mp_real_small_sc_add_2(ulong * S1, ulong * S0, ulong * G1, ulong * G0,
    ulong sa1, ulong sa0, ulong ga1, ulong ga0,
    ulong sb1, ulong sb0, ulong gb1, ulong gb0)
{
    ulong p1, p0, q1, q0, t1, t0, u1, u0, w1, w0;

    _mp_real_mulhi_2x2_sloppy(&p1, &p0, sa1, sa0, gb1, gb0);
    _mp_real_mulhi_2x2_sloppy(&q1, &q0, ga1, ga0, sb1, sb0);
    _mp_real_mulhi_2x2_sloppy(&t1, &t0, ga1, ga0, gb1, gb0);
    _mp_real_mulhi_2x2_sloppy(&u1, &u0, sa1, sa0, sb1, sb0);

    add_ssaaaa(w1, w0, sa1, sa0, sb1, sb0);
    sub_ddmmss(w1, w0, w1, w0, p1, p0);
    sub_ddmmss(*S1, *S0, w1, w0, q1, q0);

    add_ssaaaa(w1, w0, ga1, ga0, gb1, gb0);
    add_ssaaaa(w1, w0, w1, w0, u1, u0);
    sub_ddmmss(*G1, *G0, w1, w0, t1, t0);
}

FLINT_FORCE_INLINE void
_mp_real_small_sc_add_3(ulong * S2, ulong * S1, ulong * S0,
    ulong * G2, ulong * G1, ulong * G0,
    ulong sa2, ulong sa1, ulong sa0, ulong ga2, ulong ga1, ulong ga0,
    ulong sb2, ulong sb1, ulong sb0, ulong gb2, ulong gb1, ulong gb0)
{
    ulong p2, p1, p0, q2, q1, q0, t2, t1, t0, u2, u1, u0, w2, w1, w0;

    _mp_real_small_mulhi_3x3(&p2, &p1, &p0, sa2, sa1, sa0, gb2, gb1, gb0);
    _mp_real_small_mulhi_3x3(&q2, &q1, &q0, ga2, ga1, ga0, sb2, sb1, sb0);
    _mp_real_small_mulhi_3x3(&t2, &t1, &t0, ga2, ga1, ga0, gb2, gb1, gb0);
    _mp_real_small_mulhi_3x3(&u2, &u1, &u0, sa2, sa1, sa0, sb2, sb1, sb0);

    add_sssaaaaaa(w2, w1, w0, sa2, sa1, sa0, sb2, sb1, sb0);
    sub_dddmmmsss(w2, w1, w0, w2, w1, w0, p2, p1, p0);
    sub_dddmmmsss(*S2, *S1, *S0, w2, w1, w0, q2, q1, q0);

    add_sssaaaaaa(w2, w1, w0, ga2, ga1, ga0, gb2, gb1, gb0);
    add_sssaaaaaa(w2, w1, w0, w2, w1, w0, u2, u1, u0);
    sub_dddmmmsss(*G2, *G1, *G0, w2, w1, w0, t2, t1, t0);
}

/* angle addition on (sin, 1 - cos) pairs at one limb */
FLINT_FORCE_INLINE void
_mp_real_small_sc_add_1(ulong * S, ulong * G, ulong sa, ulong ga,
    ulong sb, ulong gb)
{
    ulong p = n_mulhi(sa, gb), q = n_mulhi(ga, sb);
    ulong t = n_mulhi(ga, gb), u = n_mulhi(sa, sb);

    *S = sa + sb - p - q;
    *G = ga + gb + u - t;
}

/* one fraction limb: the bounds of the two-limb kernel (the same
   operations at one limb, every product within one ulp) */
#define MP_REAL_SMALL_SIN_COS_ERR_1 32

/*
    (s, g) = (sin v, 1 - cos v) at one fraction limb for v in [0, 1):
    v = j1/2^7 + j2/2^14 + r, the table pairs combined by angle addition,
    then with sin r = r - r^3/6 and 1 - cos r = q/2 - q^2/24, q = r^2
    (the next terms below 2^-76 resp. 2^-93).
*/
FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_g_1(ulong * s, ulong * g, ulong v)
{
    ulong j1, j2, R, q, sr, gr, ts, tg;
    const ulong * T1;
    const ulong * T2;

    j1 = v >> 57;
    j2 = (v >> 50) & 127;
    R = v & ((UWORD(1) << 50) - 1);

    T1 = _mp_real_small_sin_cos_s1 + 6 * j1;
    T2 = _mp_real_small_sin_cos_s2 + 6 * j2;
    _mp_real_small_sc_add_1(&ts, &tg, T1[2], T1[5], T2[2], T2[5]);

    q = n_mulhi(R, R);
    sr = R - n_mulhi(n_mulhi(q, R), MP_REAL_SMALL_F6_2);
    gr = (q >> 1) - n_mulhi(n_mulhi(q, q), MP_REAL_SMALL_F24_2);

    _mp_real_small_sc_add_1(s, g, ts, tg, sr, gr);
}

/*
    (s, g) = (sin v, 1 - cos v) at two fraction limbs for v = (v1, v0)
    in [0, 1):

        v = j1/2^7 + j2/2^14 + r,  r < 2^-14,

    the table pairs for j1/2^7 and j2/2^14 combined by angle addition,
    then with the polynomials (q = r^2 < 2^-28)

        sin r = r - r^3 (1/6 - q (1/120 - q/5040))       (next term below 2^-144)
        1 - cos r = q/2 - q^2 (1/24 - q (1/720 - q/8!))   (next below 2^-161)

    at one limb for the innermost terms (on q 2^28), then at two.
*/
FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_g_2(ulong * s1, ulong * s0, ulong * g1, ulong * g0,
    ulong v1, ulong v0)
{
    ulong j1, j2, r1, r0, q1, q0, rho, t, p1, p0, w1, w0, sg1, sg0, gm1, gm0;
    ulong ts1, ts0, tg1, tg0;
    const ulong * T1;
    const ulong * T2;

    j1 = v1 >> 57;
    j2 = (v1 >> 50) & 127;
    r1 = v1 & ((UWORD(1) << 50) - 1);
    r0 = v0;

    /* the tables combined (independent of the polynomials) */
    T1 = _mp_real_small_sin_cos_s1 + 6 * j1;
    T2 = _mp_real_small_sin_cos_s2 + 6 * j2;
    _mp_real_small_sc_add_2(&ts1, &ts0, &tg1, &tg0,
        T1[2], T1[1], T1[5], T1[4], T2[2], T2[1], T2[5], T2[4]);

    _mp_real_sqrhi_2x2(&q1, &q0, r1, r0);
    rho = (q1 << 28) | (q0 >> 36);              /* floor(q 2^92) */

    /* sin r */
    t = MP_REAL_SMALL_F120_2 - (n_mulhi(rho, MP_REAL_SMALL_F5040_2) >> 28);
    _mp_real_mulhi_2x1(&p1, &p0, q1, q0, t);
    sub_ddmmss(p1, p0, MP_REAL_SMALL_F6_2, MP_REAL_SMALL_F6_1, p1, p0);
    _mp_real_mulhi_2x2_sloppy(&w1, &w0, q1, q0, r1, r0);
    _mp_real_mulhi_2x2_sloppy(&w1, &w0, w1, w0, p1, p0);
    sub_ddmmss(sg1, sg0, r1, r0, w1, w0);

    /* 1 - cos r */
    t = MP_REAL_SMALL_F720_2 - (n_mulhi(rho, MP_REAL_SMALL_F8_2) >> 28);
    _mp_real_mulhi_2x1(&p1, &p0, q1, q0, t);
    sub_ddmmss(p1, p0, MP_REAL_SMALL_F24_2, MP_REAL_SMALL_F24_1, p1, p0);
    _mp_real_sqrhi_2x2(&w1, &w0, q1, q0);
    _mp_real_mulhi_2x2_sloppy(&w1, &w0, w1, w0, p1, p0);
    gm1 = q1 >> 1;
    gm0 = (q0 >> 1) | (q1 << 63);
    sub_ddmmss(gm1, gm0, gm1, gm0, w1, w0);

    _mp_real_small_sc_add_2(s1, s0, g1, g0, ts1, ts0, tg1, tg0, sg1, sg0, gm1, gm0);
}

/*
    The same at three fraction limbs, with the longer polynomials

        sin r = r - r^3 (1/6 - q (1/120 - q (1/5040 - q (1/9! - q/11!))))
        1 - cos r = q/2 - q^2 (1/24 - q (1/720 - q (1/8! - q/10!)))

    (next terms below 2^-214 resp. 2^-196), at two limbs (on the top two
    limbs of q) except for the outermost steps.
*/
FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_g_3(ulong * s2, ulong * s1, ulong * s0,
    ulong * g2, ulong * g1, ulong * g0, ulong v2, ulong v1, ulong v0)
{
    ulong j1, j2, r2, r1, r0, q2, q1, q0, x1, x0, p2, p1, p0, w2, w1, w0;
    ulong sg2, sg1, sg0, gm2, gm1, gm0, ts2, ts1, ts0, tg2, tg1, tg0;
    const ulong * T1;
    const ulong * T2;

    j1 = v2 >> 57;
    j2 = (v2 >> 50) & 127;
    r2 = v2 & ((UWORD(1) << 50) - 1);
    r1 = v1;
    r0 = v0;

    T1 = _mp_real_small_sin_cos_s1 + 6 * j1;
    T2 = _mp_real_small_sin_cos_s2 + 6 * j2;
    _mp_real_small_sc_add_3(&ts2, &ts1, &ts0, &tg2, &tg1, &tg0,
        T1[2], T1[1], T1[0], T1[5], T1[4], T1[3],
        T2[2], T2[1], T2[0], T2[5], T2[4], T2[3]);

    _mp_real_small_sqrhi_3x3(&q2, &q1, &q0, r2, r1, r0);

    /* sin r */
    _mp_real_mulhi_2x1(&x1, &x0, q2, q1, MP_REAL_SMALL_F11_2);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_F9_2, MP_REAL_SMALL_F9_1, x1, x0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, q2, q1, x1, x0);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_F5040_2, MP_REAL_SMALL_F5040_1, x1, x0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, q2, q1, x1, x0);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_F120_2, MP_REAL_SMALL_F120_1, x1, x0);
    _mp_real_mulhi_3x2(&p2, &p1, &p0, q2, q1, q0, x1, x0);
    sub_dddmmmsss(p2, p1, p0, MP_REAL_SMALL_F6_2, MP_REAL_SMALL_F6_1, MP_REAL_SMALL_F6_0, p2, p1, p0);
    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, q2, q1, q0, r2, r1, r0);
    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, w2, w1, w0, p2, p1, p0);
    sub_dddmmmsss(sg2, sg1, sg0, r2, r1, r0, w2, w1, w0);

    /* 1 - cos r */
    _mp_real_mulhi_2x1(&x1, &x0, q2, q1, MP_REAL_SMALL_F10_2);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_F8_2, MP_REAL_SMALL_F8_1, x1, x0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, q2, q1, x1, x0);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_F720_2, MP_REAL_SMALL_F720_1, x1, x0);
    _mp_real_mulhi_3x2(&p2, &p1, &p0, q2, q1, q0, x1, x0);
    sub_dddmmmsss(p2, p1, p0, MP_REAL_SMALL_F24_2, MP_REAL_SMALL_F24_1, MP_REAL_SMALL_F24_0, p2, p1, p0);
    _mp_real_small_sqrhi_3x3(&w2, &w1, &w0, q2, q1, q0);
    _mp_real_small_mulhi_3x3(&w2, &w1, &w0, w2, w1, w0, p2, p1, p0);
    gm2 = q2 >> 1;
    gm1 = (q1 >> 1) | (q2 << 63);
    gm0 = (q0 >> 1) | (q1 << 63);
    sub_dddmmmsss(gm2, gm1, gm0, gm2, gm1, gm0, w2, w1, w0);

    _mp_real_small_sc_add_3(s2, s1, s0, g2, g1, g0, ts2, ts1, ts0, tg2, tg1, tg0,
        sg2, sg1, sg0, gm2, gm1, gm0);
}

/* array wrappers: (ys, n + 1), (yc, n + 1) = sin v, cos v, v at n = 1,
   2 resp. 3 fraction limbs */
FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_1_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
{
    ulong g;
    _mp_real_small_sin_cos_g_1(ys, &g, v[0]);
    ys[1] = 0;
    yc[0] = -g;
    yc[1] = (g == 0);
    *err = MP_REAL_SMALL_SIN_COS_ERR_1;
}

FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_2_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
{
    ulong g1, g0;
    _mp_real_small_sin_cos_g_2(ys + 1, ys, &g1, &g0, v[1], v[0]);
    ys[2] = 0;
    sub_dddmmmsss(yc[2], yc[1], yc[0], UWORD(1), UWORD(0), UWORD(0), UWORD(0), g1, g0);
    *err = MP_REAL_SMALL_SIN_COS_ERR_2;
}

FLINT_FORCE_INLINE void
_mp_real_small_sin_cos_3_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
{
    ulong g2, g1, g0;
    _mp_real_small_sin_cos_g_3(ys + 2, ys + 1, ys, &g2, &g1, &g0, v[2], v[1], v[0]);
    ys[3] = 0;
    sub_ddddmmmmssss(yc[3], yc[2], yc[1], yc[0], UWORD(1), UWORD(0), UWORD(0), UWORD(0),
        UWORD(0), g2, g1, g0);
    *err = MP_REAL_SMALL_SIN_COS_ERR_3;
}

/* atan ************************************************************************/

/* 1/(2k + 1), top limbs first */
#define MP_REAL_SMALL_I3 UWORD(0x5555555555555555)      /* every limb */
#define MP_REAL_SMALL_I5 UWORD(0x3333333333333333)      /* every limb */
#define MP_REAL_SMALL_I7_2 UWORD(0x2492492492492492)
#define MP_REAL_SMALL_I7_1 UWORD(0x4924924924924924)
#define MP_REAL_SMALL_I7_0 UWORD(0x9249249249249249)
#define MP_REAL_SMALL_I9_2 UWORD(0x1c71c71c71c71c71)
#define MP_REAL_SMALL_I9_1 UWORD(0xc71c71c71c71c71c)
#define MP_REAL_SMALL_I11_2 UWORD(0x1745d1745d1745d1)
#define MP_REAL_SMALL_I11_1 UWORD(0x745d1745d1745d17)
#define MP_REAL_SMALL_I13_2 UWORD(0x13b13b13b13b13b1)
#define MP_REAL_SMALL_I13_1 UWORD(0x3b13b13b13b13b13)
#define MP_REAL_SMALL_I15 UWORD(0x1111111111111111)
#define MP_REAL_SMALL_I17 UWORD(0x0f0f0f0f0f0f0f0f)
#define MP_REAL_SMALL_I19 UWORD(0x0d79435e50d79435)

/* error bounds in ulps of B^-n, the same for atan(v) and atan(Y / X):
   at two limbs w within one ulp (the floor of an exact quotient, or for
   atan(v) of the quotient by a truncated D, which only raises it by
   w 2^-128), w^3 P within 4, A within 1, so within 6; at three limbs w
   within 9 (_mp_real_small_divw_3), w^3 P within 6.7, A within 1 and the
   dropped term w^21/21 below 1/2, so within 17.2 */
#define MP_REAL_SMALL_ATAN_ERR_2 16
#define MP_REAL_SMALL_ATAN_ERR_3 32

/*
    The arctangent kernels reduce by c = j/2^9,

        atan(t) = A[j] + atan(w),  w = (t - c) / (1 + t c) in [0, 2^-9),

    with w by one division in registers: for t = v at n fraction limbs,
    w = (v - c) / (1 + v c); for a ratio t = Y / X, w = (512 Y - j X) /
    (512 X + j Y) with exact numerator and denominator (and j = 512 for
    Y = X). Then atan(w) = w - w^3 P(z), z = w^2, by the odd series (the
    next term below 2^-138 resp. 2^-193 for n = 2, 3), the innermost terms
    at one limb (on z 2^18), then at two limbs while the powers of z allow,
    then at three.
*/

/* (y1, y0) = A[j] + atan(w) at two fraction limbs, w = (w1, w0) < 2^-9 */
FLINT_FORCE_INLINE void
_mp_real_small_atan_w_2(ulong * y1, ulong * y0, ulong w1, ulong w0, ulong j)
{
    ulong z1, z0, zeta, t, r1, r0;
    const ulong * A;

    /* z = w^2 < 2^-18; S = 1/9 - z/11 + z^2/13 at one limb */
    _mp_real_sqrhi_2x2(&z1, &z0, w1, w0);
    zeta = (z1 << 18) | (z0 >> 46);         /* floor(z 2^82) */
    t = MP_REAL_SMALL_I11_2 - (n_mulhi(zeta, MP_REAL_SMALL_I13_2) >> 18);
    t = MP_REAL_SMALL_I9_2 - (n_mulhi(zeta, t) >> 18);

    /* R = 1/7 - z S, Q = 1/5 - z R, P = 1/3 - z Q */
    _mp_real_mulhi_2x1(&r1, &r0, z1, z0, t);
    sub_ddmmss(r1, r0, MP_REAL_SMALL_I7_2, MP_REAL_SMALL_I7_1, r1, r0);
    _mp_real_mulhi_2x2_sloppy(&r1, &r0, z1, z0, r1, r0);
    sub_ddmmss(r1, r0, MP_REAL_SMALL_I5, MP_REAL_SMALL_I5, r1, r0);
    _mp_real_mulhi_2x2_sloppy(&r1, &r0, z1, z0, r1, r0);
    sub_ddmmss(r1, r0, MP_REAL_SMALL_I3, MP_REAL_SMALL_I3, r1, r0);

    /* atan(w) = w - w^3 P */
    _mp_real_mulhi_2x2_sloppy(&z1, &z0, z1, z0, w1, w0);
    _mp_real_mulhi_2x2_sloppy(&z1, &z0, z1, z0, r1, r0);
    sub_ddmmss(w1, w0, w1, w0, z1, z0);

    A = _mp_real_small_atan_a + 3 * j;
    add_ssaaaa(*y1, *y0, A[2], A[1], w1, w0);
}

/* (y2, y1, y0) = A[j] + atan(w) at three fraction limbs, w < 2^-9 */
FLINT_FORCE_INLINE void
_mp_real_small_atan_w_3(ulong * y2, ulong * y1, ulong * y0, ulong w2, ulong w1, ulong w0, ulong j)
{
    ulong z2, z1, z0, zeta, t, x1, x0, r2, r1, r0;
    const ulong * A;

    /* w from _mp_real_small_divw_3 may exceed the exact value (below
       2^-9) by a few ulps: clamp it below 2^-9, which only moves it
       towards the exact value */
    if (FLINT_UNLIKELY(w2 >= (UWORD(1) << 55)))
    {
        w2 = (UWORD(1) << 55) - 1;
        w1 = w0 = UWORD_MAX;
    }

    /* z = w^2 < 2^-18; 1/15 - z (1/17 - z/19) at one limb */
    _mp_real_small_sqrhi_3x3(&z2, &z1, &z0, w2, w1, w0);
    zeta = (z2 << 18) | (z1 >> 46);         /* floor(z 2^82) */
    t = MP_REAL_SMALL_I17 - (n_mulhi(zeta, MP_REAL_SMALL_I19) >> 18);
    t = MP_REAL_SMALL_I15 - (n_mulhi(zeta, t) >> 18);

    /* 1/13 - z (..), 1/11 - z (..), 1/9 - z (..) at two limbs */
    _mp_real_mulhi_2x1(&x1, &x0, z2, z1, t);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_I13_2, MP_REAL_SMALL_I13_1, x1, x0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, z2, z1, x1, x0);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_I11_2, MP_REAL_SMALL_I11_1, x1, x0);
    _mp_real_mulhi_2x2_sloppy(&x1, &x0, z2, z1, x1, x0);
    sub_ddmmss(x1, x0, MP_REAL_SMALL_I9_2, MP_REAL_SMALL_I9_1, x1, x0);

    /* 1/7 - z (..), 1/5 - z (..), 1/3 - z (..) at three limbs */
    _mp_real_mulhi_3x2(&r2, &r1, &r0, z2, z1, z0, x1, x0);
    sub_dddmmmsss(r2, r1, r0, MP_REAL_SMALL_I7_2, MP_REAL_SMALL_I7_1, MP_REAL_SMALL_I7_0, r2, r1, r0);
    _mp_real_small_mulhi_3x3(&r2, &r1, &r0, z2, z1, z0, r2, r1, r0);
    sub_dddmmmsss(r2, r1, r0, MP_REAL_SMALL_I5, MP_REAL_SMALL_I5, MP_REAL_SMALL_I5, r2, r1, r0);
    _mp_real_small_mulhi_3x3(&r2, &r1, &r0, z2, z1, z0, r2, r1, r0);
    sub_dddmmmsss(r2, r1, r0, MP_REAL_SMALL_I3, MP_REAL_SMALL_I3, MP_REAL_SMALL_I3, r2, r1, r0);

    /* atan(w) = w - w^3 P */
    _mp_real_small_mulhi_3x3(&z2, &z1, &z0, z2, z1, z0, w2, w1, w0);
    _mp_real_small_mulhi_3x3(&z2, &z1, &z0, z2, z1, z0, r2, r1, r0);
    sub_dddmmmsss(w2, w1, w0, w2, w1, w0, z2, z1, z0);

    A = _mp_real_small_atan_a + 3 * j;
    add_sssaaaaaa(*y2, *y1, *y0, A[2], A[1], A[0], w2, w1, w0);
}

/*
    (w2, w1, w0) = N / D at three fraction limbs for a fraction N =
    (n2, n1, n0) and D = dI + (d2, d1, d0), 1 <= dI < 2^10, with N / D <
    2^-9: w' = (q1, q0, 0) the two-limb quotient of the top limbs (within
    2^-127 of N / D), then w = w' + e / D for the residual e = N - w' D
    (within 3 ulps above, |e| < 2^-117) times a one-limb inverse of D
    (relative error below 2^-63), within 9 ulps in all.
*/
FLINT_FORCE_INLINE void
_mp_real_small_divw_3(ulong * w2, ulong * w1, ulong * w0,
    ulong n2, ulong n1, ulong n0, ulong dI, ulong d2, ulong d1, ulong d0)
{
    ulong q[2], e2, e1, e0, s2, s1, s0, h, l, p2, p1, dt, inv, rem;
    unsigned int sh, k;
    int neg;

    _mp_real_divq_4_3z(q, n2, n1, dI, d2, d1);

    /* e = N - w' dI - w' d, the first product exact (below N) */
    umul_ppmm(p2, p1, q[0], dI);
    umul_ppmm(h, l, q[1], dI);
    add_ssaaaa(h, p2, h, p2, UWORD(0), l);
    (void) h;
    _mp_real_mulhi_3x2(&s2, &s1, &s0, d2, d1, d0, q[1], q[0]);
    sub_dddmmmsss(e2, e1, e0, n2, n1, n0, p2, p1, UWORD(0));
    sub_dddmmmsss(e2, e1, e0, e2, e1, e0, s2, s1, s0);
    neg = (slong) e2 < 0;
    if (neg)
        sub_dddmmmsss(e2, e1, e0, UWORD(0), UWORD(0), UWORD(0), e2, e1, e0);
    (void) e2;

    /* inv = floor((2^127 - 1) / Dt) ~ 2^(127 - sh) / D, Dt the top 64
       bits of D 2^sh */
    sh = flint_clz(dI);
    dt = (dI << sh) | (d2 >> (FLINT_BITS - sh));
    udiv_qrnnd(inv, rem, (UWORD(1) << 63) - 1, UWORD_MAX, dt);
    (void) rem;

    /* e / D = (e inv / B) 2^(sh - 63) */
    _mp_real_mulhi_2x1(&s1, &s0, e1, e0, inv);
    k = (FLINT_BITS - 1) - sh;
    if (k != 0)
    {
        s0 = (s0 >> k) | (s1 << (FLINT_BITS - k));
        s1 = s1 >> k;
    }

    if (neg)
        sub_dddmmmsss(*w2, *w1, *w0, q[1], q[0], UWORD(0), UWORD(0), s1, s0);
    else
        add_sssaaaaaa(*w2, *w1, *w0, q[1], q[0], UWORD(0), UWORD(0), s1, s0);
}

/*
    (w2, w1, w0) = floor(N B^3 / D), N / D at three fraction limbs below
    by under one ulp, for fractions N < D, D = (d2, d1, d0) normalized
    (d2 >= B/2): schoolbook with 3/2 divisions by the top two limbs of D,
    each quotient limb at most one too large and corrected, as in
    _mp_real_divq_4_3z.
*/
#define MP_REAL_SMALL_DIVQ_3_ERR 1

FLINT_FORCE_INLINE void
_mp_real_small_divq_3(ulong * w2, ulong * w1, ulong * w0,
    ulong n2, ulong n1, ulong n0, ulong d2, ulong d1, ulong d0)
{
    ulong u2 = n2, u1 = n1, u0 = n0, qq, r1, r0, p1, p0, dinv, q[3];
    int j;

    dinv = _MP_REAL_PREINV1(d2, d1);
    (void) dinv;

    /* the partial remainder (u2, u1, u0) < D, followed by zero limbs */
    for (j = 2; j >= 0; j--)
    {
        if (FLINT_UNLIKELY(u2 == d2 && u1 == d1))
        {
            /* (the difference goes through a temporary: sub_dddmmmsss
               may write u0 before it reads its other inputs) */
            ulong t = d0 - u0;
            qq = ~UWORD(0);
            sub_dddmmmsss(u2, u1, u0, d2, d1, d0, UWORD(0), t, UWORD(0));
        }
        else
        {
            int neg;

            _MP_REAL_UDIV_QR_3BY2(qq, r1, r0, u2, u1, u0, d2, d1, dinv);
            umul_ppmm(p1, p0, qq, d0);
            neg = (r1 == 0) & ((r0 < p1) | ((r0 == p1) & (p0 != 0)));
            sub_dddmmmsss(u2, u1, u0, r1, r0, UWORD(0), UWORD(0), p1, p0);
            if (FLINT_UNLIKELY(neg))
            {
                qq--;
                add_sssaaaaaa(u2, u1, u0, u2, u1, u0, d2, d1, d0);
            }
        }

        q[j] = qq;
    }

    *w2 = q[2];
    *w1 = q[1];
    *w0 = q[0];
}

/* one fraction limb: w within one ulp either way (the floor of a
   quotient by a denominator truncated by under 2^-62), w^3 P below by
   under 1.5, A below by under 1, the dropped w^9/9 below 2^-84 */
#define MP_REAL_SMALL_ATAN_ERR_1 8

/*
    atan(v) at one fraction limb for v in [0, 1): with j = floor(512 v),
    w = (v - j/512) / (1 + v j/512) < 2^-9 by one hardware division, and
    atan(v) = A[j] + w - w^3 (1/3 - z (1/5 - z/7)), z = w^2.
*/
FLINT_FORCE_INLINE ulong
_mp_real_small_atan_1(ulong v)
{
    ulong j, n, h, l, d, w, rem, z, p;

    j = v >> 55;
    n = v & ((UWORD(1) << 55) - 1);

    /* D = 1 + d, d = v j / 2^9 truncated; Dn = D 2^63 */
    umul_ppmm(h, l, v, j);
    d = (h << 55) | (l >> 9);
    udiv_qrnnd(w, rem, n >> 1, n << 63, (UWORD(1) << 63) | (d >> 1));
    (void) rem;

    z = n_mulhi(w, w);
    p = MP_REAL_SMALL_I5 - n_mulhi(z, MP_REAL_SMALL_I7_2);
    p = MP_REAL_SMALL_I3 - n_mulhi(z, p);
    return _mp_real_small_atan_a[3 * j + 2] + w - n_mulhi(n_mulhi(z, w), p);
}

/* (y1, y0) = atan(v) at two fraction limbs for v = (v1, v0) in [0, 1) */
FLINT_FORCE_INLINE void
_mp_real_small_atan_2(ulong * y1, ulong * y0, ulong v1, ulong v0)
{
    ulong j, n1, n0, p2, p1, p0, h, l, d1, d0, q[2];

    j = v1 >> 55;
    n1 = v1 & ((UWORD(1) << 55) - 1);
    n0 = v0;

    /* D = 1 + v j / 2^9, the fraction truncated to two limbs */
    umul_ppmm(p1, p0, v0, j);
    umul_ppmm(h, l, v1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    d1 = (p2 << 55) | (p1 >> 9);
    d0 = (p1 << 55) | (p0 >> 9);

    _mp_real_divq_4_3z(q, n1, n0, UWORD(1), d1, d0);
    _mp_real_small_atan_w_2(y1, y0, q[1], q[0], j);
}

/* (y2, y1, y0) = atan(v) at three fraction limbs */
FLINT_FORCE_INLINE void
_mp_real_small_atan_3(ulong * y2, ulong * y1, ulong * y0, ulong v2, ulong v1, ulong v0)
{
    ulong j, n2, p3, p2, p1, p0, h, l, d2, d1, d0, w2, w1, w0;

    j = v2 >> 55;
    n2 = v2 & ((UWORD(1) << 55) - 1);

    /* D = 1 + v j / 2^9, the fraction truncated to three limbs */
    umul_ppmm(p1, p0, v0, j);
    umul_ppmm(h, l, v1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, v2, j);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    d2 = (p3 << 55) | (p2 >> 9);
    d1 = (p2 << 55) | (p1 >> 9);
    d0 = (p1 << 55) | (p0 >> 9);

    _mp_real_small_divw_3(&w2, &w1, &w0, n2, v1, v0, UWORD(1), d2, d1, d0);
    _mp_real_small_atan_w_3(y2, y1, y0, w2, w1, w0, j);
}

/* j = floor(512 Y / X), or one off, by a hardware division of the top
   limbs (512 yt / xt within 2^-53 of 512 Y / X for X >= 1/2) */
FLINT_FORCE_INLINE ulong
_mp_real_small_atan_ratio_j(ulong yt, ulong xt)
{
    ulong j, r;
    udiv_qrnnd(j, r, yt >> 55, yt << 9, xt);
    (void) r;
    return j;
}

/*
    (r1, r0) = atan(Y / X) at two fraction limbs for fractions X = (x1, x0)
    >= 1/2 and Y = (y1, y0) <= X, by one division: with j = floor(512 Y /
    X) (from the top limbs, corrected exactly),
    w = (512 Y - j X) / (512 X + j Y). (One hardware division more than
    the kernel for atan(v), for j, against two resp. three for a
    reciprocal followed by atan(v).)
*/
FLINT_FORCE_INLINE void
_mp_real_small_atan_ratio_2(ulong * r1, ulong * r0, ulong y1, ulong y0, ulong x1, ulong x0)
{
    ulong j, n2, n1, n0, p2, p1, p0, h, l, dI, d1, d0, q[2];

    j = _mp_real_small_atan_ratio_j(y1, x1);

    /* N = 512 Y - j X, exactly, with an integral limb */
    n2 = y1 >> 55;
    n1 = (y1 << 9) | (y0 >> 55);
    n0 = y0 << 9;
    umul_ppmm(p1, p0, x0, j);
    umul_ppmm(h, l, x1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    sub_dddmmmsss(n2, n1, n0, n2, n1, n0, p2, p1, p0);

    /* now 0 <= N < X */
    if ((slong) n2 < 0)
    {
        j--;
        add_sssaaaaaa(n2, n1, n0, n2, n1, n0, UWORD(0), x1, x0);
    }
    else if (n2 != 0 || n1 > x1 || (n1 == x1 && n0 >= x0))
    {
        j++;
        sub_dddmmmsss(n2, n1, n0, n2, n1, n0, UWORD(0), x1, x0);
    }
    FLINT_ASSERT(n2 == 0 && j <= 512);

    /* D = 512 X + j Y in [256, 1024) */
    dI = x1 >> 55;
    d1 = (x1 << 9) | (x0 >> 55);
    d0 = x0 << 9;
    umul_ppmm(p1, p0, y0, j);
    umul_ppmm(h, l, y1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    add_sssaaaaaa(dI, d1, d0, dI, d1, d0, p2, p1, p0);

    _mp_real_divq_4_3z(q, n1, n0, dI, d1, d0);
    _mp_real_small_atan_w_2(r1, r0, q[1], q[0], j);
}

/* the same at three fraction limbs */
FLINT_FORCE_INLINE void
_mp_real_small_atan_ratio_3(ulong * r2, ulong * r1, ulong * r0,
    ulong y2, ulong y1, ulong y0, ulong x2, ulong x1, ulong x0)
{
    ulong j, n3, n2, n1, n0, p3, p2, p1, p0, h, l, dI, d2, d1, d0, w2, w1, w0;

    j = _mp_real_small_atan_ratio_j(y2, x2);

    n3 = y2 >> 55;
    n2 = (y2 << 9) | (y1 >> 55);
    n1 = (y1 << 9) | (y0 >> 55);
    n0 = y0 << 9;
    umul_ppmm(p1, p0, x0, j);
    umul_ppmm(h, l, x1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, x2, j);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    sub_ddddmmmmssss(n3, n2, n1, n0, n3, n2, n1, n0, p3, p2, p1, p0);

    if ((slong) n3 < 0)
    {
        j--;
        add_ssssaaaaaaaa(n3, n2, n1, n0, n3, n2, n1, n0, UWORD(0), x2, x1, x0);
    }
    else if (n3 != 0 || n2 > x2 || (n2 == x2 && (n1 > x1 || (n1 == x1 && n0 >= x0))))
    {
        j++;
        sub_ddddmmmmssss(n3, n2, n1, n0, n3, n2, n1, n0, UWORD(0), x2, x1, x0);
    }
    FLINT_ASSERT(n3 == 0 && j <= 512);

    dI = x2 >> 55;
    d2 = (x2 << 9) | (x1 >> 55);
    d1 = (x1 << 9) | (x0 >> 55);
    d0 = x0 << 9;
    umul_ppmm(p1, p0, y0, j);
    umul_ppmm(h, l, y1, j);
    add_ssaaaa(p2, p1, h, l, UWORD(0), p1);
    umul_ppmm(h, l, y2, j);
    add_ssaaaa(p3, p2, h, l, UWORD(0), p2);
    add_ssssaaaaaaaa(dI, d2, d1, d0, dI, d2, d1, d0, p3, p2, p1, p0);

    _mp_real_small_divw_3(&w2, &w1, &w0, n2, n1, n0, dI, d2, d1, d0);
    _mp_real_small_atan_w_3(r2, r1, r0, w2, w1, w0, j);
}

/* array wrappers: (y, n) = atan(v), v at n = 1, 2 resp. 3 fraction
   limbs */
FLINT_FORCE_INLINE void
_mp_real_small_atan_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    y[0] = _mp_real_small_atan_1(v[0]);
    *err = MP_REAL_SMALL_ATAN_ERR_1;
}

FLINT_FORCE_INLINE void
_mp_real_small_atan_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_atan_2(y + 1, y, v[1], v[0]);
    *err = MP_REAL_SMALL_ATAN_ERR_2;
}

FLINT_FORCE_INLINE void
_mp_real_small_atan_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
{
    _mp_real_small_atan_3(y + 2, y + 1, y, v[2], v[1], v[0]);
    *err = MP_REAL_SMALL_ATAN_ERR_3;
}

#endif /* FLINT_BITS == 64 */

#endif

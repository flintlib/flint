/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

/* Fast implementation of the Sorenson-Webster certified Rabin-Miller
   primality test for integers up to about 81 bits. More precisely, the
   test is valid for n < SWbound = 3317044064679887385961981. */

static const ulong SWbound[2] = { UWORD(0x51adc5b22410a5fd), UWORD(0x2be69) };


/* GCC 11/Zen 3: the following functions are faster implemented with
   umul_ppmm/add_ssaaaa inline assembly than with __uint128_t. */

FLINT_FORCE_INLINE
void sqrlo_3x2x2(ulong * r2, ulong * r1, ulong * r0, ulong a1, ulong a0)
{
    mp_limb_t t1, t2, s0, s1, s2;
    umul_ppmm(s1, s0, a0, a0);
    s2 = a1 * a1;
    /* Adding twice is faster than shifting. When we know that 2*a1 does
       not overflow, we could just do a0*(2*a1), but this is also
       (surprisingly) slower.  */
    umul_ppmm(t2, t1, a0, a1);
    add_ssaaaa(s2, s1, s2, s1, t2, t1);
    add_ssaaaa(s2, s1, s2, s1, t2, t1);
    *r0 = s0;
    *r1 = s1;
    *r2 = s2;
}

FLINT_FORCE_INLINE
void mulhi_2x2_sloppy(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b1, ulong b0)
{
    ulong s2, u2, v3, v2, r3, r2;
    s2 = n_mulhi(a0, b1);
    u2 = n_mulhi(a1, b0);
    add_ssaaaa(r3, r2, 0, s2, 0, u2);
    umul_ppmm(v3, v2, a1, b1);
    add_ssaaaa(r3, r2, r3, r2, v3, v2);
    *r0 = r2;
    *r1 = r3;
}

FLINT_FORCE_INLINE
void mullo_2x1(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b0)
{
    ulong t0, t1;
    umul_ppmm(t1, t0, a0, b0);
    t1 += a1 * b0;
    *r0 = t0;
    *r1 = t1;
}

FLINT_FORCE_INLINE
void mullo_2x2(ulong * r1, ulong * r0, ulong a1, ulong a0, ulong b1, ulong b0)
{
    FLINT_MPN_MULLOW_2X2(*r1, *r0, a1, a0, b1, b0);
}

/*
n_ll_small: "1.5-limb" arithmetic mod m where B < m < B^(3/2), B = 2^64.
We use the approximate inverse minv = floor(B^3 / m). Given 0 <= x < B^3,

    q' = mulhi_2x2_sloppy(floor(x / B), minv).

is an approximation of the true floor quotient q = floor(x / m). We have
0 <= q' < B^2, and

    max(0, q-5) <= q' <= q.

An approximate modular reduction using q' gives a noncanonical
residue satisfying

    0  <=  x - q' * m  <  6m.

This is useful for Miller-Rabin tests, e.g. in our application where we want
to compute b^e mod m with m < 2^82 and b <= 41, a valid noncanonical residue
y < 6m as input to a multiply-and-square step gives

    (y*b)^2 <= (6 * 2^82 * 41)^2 < 2^192.

*/

/* Compute minv = 2^192 / m. TODO: can we beat mpn_tdiv_qr here? */
void n_ll_small_preinv(nn_ptr minv, nn_srcptr m)
{
    ulong qq[3];
    ulong nn[4] = { 0, 0, 0, 1 } ;
    ulong rr[2];
    mpn_tdiv_qr(qq, rr, 0, nn, 4, m, 2);
    minv[0] = qq[0];
    minv[1] = qq[1];
}

/* Reduce three-limb value {x0,x1,x2} to noncanonical residue {r0,r1} */
FLINT_FORCE_INLINE
void n_ll_small_reduce3_sloppy(ulong * r1, ulong *r0,
    ulong x2, ulong x1, ulong x0,
    ulong m1, ulong m0, ulong minv1, ulong minv0)
{
    ulong q1, q0, t1, t0;
    mulhi_2x2_sloppy(&q1, &q0, x2, x1, minv1, minv0);
    mullo_2x2(&t1, &t0, q1, q0, m1, m0);
    sub_ddmmss(x1, x0, x1, x0, t1, t0);
    *r0 = x0;
    *r1 = x1;
}

/* Reduce two-limb value to a noncanonical residue.
   The result will typically be in [0, 2m). */
FLINT_FORCE_INLINE
void n_ll_small_reduce2_sloppy(ulong * r1, ulong * r0,
    ulong x1, ulong x0, ulong m1, ulong m0, ulong minv1, ulong FLINT_UNUSED(minv0))
{
    ulong q0, t1, t0;
    q0 = n_mulhi(x1, minv1);
    mullo_2x1(&t1, &t0, m1, m0, q0);
    sub_ddmmss(x1, x0, x1, x0, t1, t0);
    *r0 = x0;
    *r1 = x1;
}

/* Reduce two-limb value to [0,m). */
FLINT_FORCE_INLINE
void n_ll_small_reduce2_precise(ulong * r1, ulong * r0,
    ulong x1, ulong x0, ulong m1, ulong m0, ulong minv1, ulong minv0)
{
    n_ll_small_reduce2_sloppy(&x1, &x0, x1, x0, m1, m0, minv1, minv0);

    /* Todo: it may be possible to show that this is needed at most once. */
    while (x1 > m1 || (x1 == m1 && x0 >= m0))
        sub_ddmmss(x1, x0, x1, x0, m1, m0);
    *r0 = x0;
    *r1 = x1;
}

/* Squaring with reduction to noncanonical residue. */
FLINT_FORCE_INLINE
void n_ll_small_sqrmod_sloppy(ulong * r1, ulong * r0, ulong a1, ulong a0,
    ulong m1, ulong m0, ulong minv1, ulong minv0)
{
    ulong t2, t1, t0, q1, q0, u1, u0;

    sqrlo_3x2x2(&t2, &t1, &t0, a1, a0);
    mulhi_2x2_sloppy(&q1, &q0, t2, t1, minv1, minv0);
    mullo_2x2(&u1, &u0, q1, q0, m1, m0);
    sub_ddmmss(t1, t0, t1, t0, u1, u0);

    *r1 = t1;
    *r0 = t0;
}

/* Squaring with reduction to canonical residue. */
FLINT_FORCE_INLINE
void n_ll_small_sqrmod(ulong * r1, ulong * r0, ulong a1, ulong a0,
    ulong m1, ulong m0, ulong minv1, ulong minv0)
{
    ulong t1, t0;

    n_ll_small_sqrmod_sloppy(&t1, &t0, a1, a0, m1, m0, minv1, minv0);

    /* The loop may need several iterations, but on average it seems
       to terminate quickly enough that doing reduce_2 sloppy first
       isn't faster. */
    /* n_ll_small_reduce2_sloppy(&t1, &t0, t1, t0, m1, m0, minv1, minv0); */

    while (t1 > m1 || (t1 == m1 && t0 >= m0))
        sub_ddmmss(t1, t0, t1, t0, m1, m0);

    *r1 = t1;
    *r0 = t0;
}

/*
Simultaneously compute b1^exp, b2^exp, b3^exp (mod m) to benefit from
instruction level parallelism. On Zen 3, processing two bases simultaneously
is 1.4x faster than doing two independent exponentations; three bases is 1.8x
faster than three independent exponentations. Increasing to four bases does
not result in a further speedup.

True SIMD might be even better, at least with AVX512.
*/
void
n_ll_small_powmod_triple(nn_ptr res1, nn_ptr res2, nn_ptr res3,
            ulong b1, ulong b2, ulong b3,
            nn_srcptr exp, nn_srcptr m, nn_srcptr minv)
{
    slong bit, limbi, elimbs, ebits;
    ulong m0, m1, minv0, minv1, ei;
    ulong xr0, xr1, yr0, yr1, zr0, zr1;

    m0 = m[0];
    m1 = m[1];
    minv0 = minv[0];
    minv1 = minv[1];

    xr0 = 1;
    xr1 = 0;
    yr0 = 1;
    yr1 = 0;
    zr0 = 1;
    zr1 = 0;

    /* The exponent has 1 or 2 limbs. */
    if (exp[1] == 0)
    {
        elimbs = 1;
        ebits = FLINT_BIT_COUNT(exp[0]);
    }
    else
    {
        elimbs = 2;
        ebits = FLINT_BITS + FLINT_BIT_COUNT(exp[1]);
    }

    /* Do the first two rounds without modular reduction if b^7 fits in a limb.
       To do: could do one more round, with a two-word result. */
#define MAX_7THROOT (FLINT_BITS == 64 ? 565 : 23)
    FLINT_ASSERT(b1 <= b2);
    FLINT_ASSERT(b2 <= b3);

    if (ebits >= 3 && b3 <= MAX_7THROOT)
    {
        ulong bit2, bit1;

        bit1 = exp[(ebits - 2) / FLINT_BITS] >> ((ebits - 2) % FLINT_BITS);
        bit2 = exp[(ebits - 3) / FLINT_BITS] >> ((ebits - 3) % FLINT_BITS);

        xr0 = b1 * b1;
        yr0 = b2 * b2;
        zr0 = b3 * b3;
        if (bit1 & 1)
        {
            xr0 *= b1;
            yr0 *= b2;
            zr0 *= b3;
        }

        xr0 = xr0 * xr0;
        yr0 = yr0 * yr0;
        zr0 = zr0 * zr0;
        if (bit2 & 1)
        {
            xr0 *= b1;
            yr0 *= b2;
            zr0 *= b3;
        }

        ebits -= 3;
        elimbs = (ebits + FLINT_BITS) / FLINT_BITS;

        if (ebits >= FLINT_BITS)
            ebits -= FLINT_BITS;
    }
    else
    {
        if (ebits > FLINT_BITS)
            ebits -= FLINT_BITS;
    }

    for (limbi = elimbs - 1; limbi >= 0; limbi--)
    {
        ei = exp[limbi];

        for (bit = ebits - 1; bit >= 0; bit--)
        {
            n_ll_small_sqrmod_sloppy(&xr1, &xr0, xr1, xr0, m1, m0, minv1, minv0);
            n_ll_small_sqrmod_sloppy(&yr1, &yr0, yr1, yr0, m1, m0, minv1, minv0);
            n_ll_small_sqrmod_sloppy(&zr1, &zr0, zr1, zr0, m1, m0, minv1, minv0);

            if ((ei >> bit) & 1)
            {
                mullo_2x1(&xr1, &xr0, xr1, xr0, b1);
                mullo_2x1(&yr1, &yr0, yr1, yr0, b2);
                mullo_2x1(&zr1, &zr0, zr1, zr0, b3);
            }
        }

        ebits = FLINT_BITS;
    }

    n_ll_small_reduce2_precise(res1 + 1, res1, xr1, xr0, m1, m0, minv1, minv0);
    n_ll_small_reduce2_precise(res2 + 1, res2, yr1, yr0, m1, m0, minv1, minv0);
    n_ll_small_reduce2_precise(res3 + 1, res3, zr1, zr0, m1, m0, minv1, minv0);
}

#if FLINT_BITS == 64
#define LG_FLINT_BITS 6
#else
#define LG_FLINT_BITS 5
#endif

/* Specialization for b = 2. */
void n_ll_small_2_powmod(nn_ptr res, nn_srcptr exp, nn_srcptr m, nn_srcptr minv)
{
    slong bit, limbi, elimbs, ebits;
    ulong m0, m1, minv0, minv1, ei;
    ulong xr0, xr1;

    m0 = m[0];
    m1 = m[1];
    minv0 = minv[0];
    minv1 = minv[1];

    xr0 = 1;
    xr1 = 0;

    if (exp[1] == 0)
    {
        elimbs = 1;
        ebits = FLINT_BIT_COUNT(exp[0]);
    }
    else
    {
        elimbs = 2;
        ebits = FLINT_BITS + FLINT_BIT_COUNT(exp[1]);
    }

    if (ebits >= LG_FLINT_BITS)
    {
        ulong exp_top;
        ulong shift = ebits - LG_FLINT_BITS;

        if (shift >= FLINT_BITS)
            exp_top = exp[1] >> (shift - FLINT_BITS);
        else if (ebits >= FLINT_BITS)
            exp_top = (exp[1] << (FLINT_BITS - shift)) | (exp[0] >> shift);
        else
            exp_top = exp[0] >> shift;

        xr0 = UWORD(1) << exp_top;
        ebits -= LG_FLINT_BITS;
        elimbs = (ebits + FLINT_BITS) / FLINT_BITS;
        if (ebits >= FLINT_BITS)
            ebits -= FLINT_BITS;
    }

    for (limbi = elimbs - 1; limbi >= 0; limbi--)
    {
        ei = exp[limbi];

        for (bit = ebits - 1; bit >= 0; bit--)
        {
            n_ll_small_sqrmod_sloppy(&xr1, &xr0, xr1, xr0, m1, m0, minv1, minv0);

            /* Zen 3: branching is faster than bit-fiddling, and
               add_ssaaaa is faster than shifting. */
            if ((ei >> bit) & 1)
                add_ssaaaa(xr1, xr0, xr1, xr0, xr1, xr0);
        }

        ebits = FLINT_BITS;
    }

    n_ll_small_reduce2_precise(res + 1, res, xr1, xr0, m1, m0, minv1, minv0);
}

#if FLINT_BITS == 32

/* Todo: could implement a reasonable 64-bit primality test here. */
int
n_ll_is_prime(ulong nhi, ulong nlo)
{
    return -1;
}

#else

static
int n_ll_small_sprp_stage2(nn_srcptr y, slong s, nn_srcptr nm1, nn_srcptr n, nn_srcptr ninv)
{
    int res = 0;
    ulong y0 = y[0], y1 = y[1];

    /* y == 1 */
    if (y1 == 0 && y0 == 1)
    {
        res = 1;
    }
    else
    {
        /* y == nm1 */
        res = (y0 == nm1[0]) && (y1 == nm1[1]);

        for (s--; s > 0 && !res; s--)
        {
            /* y = y^2 mod n */
            n_ll_small_sqrmod(&y1, &y0, y1, y0, n[1], n[0], ninv[1], ninv[0]);

            res = (y0 == nm1[0]) && (y1 == nm1[1]);
        }
    }

    return res;
}

static int
n_ll_small_is_prime_mr(ulong n1, ulong n0, int base2only)
{
    ulong nm1[2], t[2], y[2];
    slong s = 0;
    slong i;
    ulong n[2] = { n0, n1 };
    ulong ninv[2];

    FLINT_ASSERT(n1 != 0);
    FLINT_ASSERT(n0 % 2 != 0);

    /* nm1 = n - 1 */
    sub_ddmmss(nm1[1], nm1[0], n1, n0, 0, 1);

    /* s = val2(nm1) */
    /* t = nm1 / 2^s */
    if (nm1[0] == 0)
    {
        s = FLINT_BITS + flint_ctz(nm1[1]);
        t[1] = 0;
        t[0] = nm1[1] >> (s - FLINT_BITS);
    }
    else
    {
        s = flint_ctz(nm1[0]);
        t[1] = nm1[1] >> s;
        t[0] = (nm1[0] >> s) | (nm1[1] << (FLINT_BITS - s));
    }

    n_ll_small_preinv(ninv, n);

    /* Base-2 test. */
    n_ll_small_2_powmod(y, t, n, ninv);
    if (!n_ll_small_sprp_stage2(y, s, nm1, n, ninv))
        return 0;

    /* Only a probable prime */
    if (base2only)
        return -1;

    /* Odd base tests. */
    static const unsigned char bases[12] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 }; 

    for (i = 0; i < 12; i += 3)
    {
        ulong y1[2], y2[2], y3[2];

        n_ll_small_powmod_triple(y1, y2, y3, bases[i], bases[i + 1], bases[i + 2], t, n, ninv);

        if (!n_ll_small_sprp_stage2(y1, s, nm1, n, ninv))
            return 0;
        if (!n_ll_small_sprp_stage2(y2, s, nm1, n, ninv))
            return 0;
        if (!n_ll_small_sprp_stage2(y3, s, nm1, n, ninv))
            return 0;
    }

    return 1;
}

/* Low and high word of 1/p mod 2^128; low and high word of (2^128-1) / p */
static const ulong ll_trial_primes[64][4] = {
    { UWORD(0xaaaaaaaaaaaaaaab), UWORD(0xaaaaaaaaaaaaaaaa), UWORD(0x5555555555555555), UWORD(0x5555555555555555) },  // 3
    { UWORD(0xcccccccccccccccd), UWORD(0xcccccccccccccccc), UWORD(0x3333333333333333), UWORD(0x3333333333333333) },  // 5
    { UWORD(0x6db6db6db6db6db7), UWORD(0xb6db6db6db6db6db), UWORD(0x4924924924924924), UWORD(0x2492492492492492) },  // 7
    { UWORD(0x2e8ba2e8ba2e8ba3), UWORD(0xa2e8ba2e8ba2e8ba), UWORD(0x745d1745d1745d17), UWORD(0x1745d1745d1745d1) },  // 11
    { UWORD(0x4ec4ec4ec4ec4ec5), UWORD(0xc4ec4ec4ec4ec4ec), UWORD(0x3b13b13b13b13b13), UWORD(0x13b13b13b13b13b1) },  // 13
    { UWORD(0xf0f0f0f0f0f0f0f1), UWORD(0xf0f0f0f0f0f0f0f0), UWORD(0xf0f0f0f0f0f0f0f), UWORD(0xf0f0f0f0f0f0f0f) },  // 17
    { UWORD(0x86bca1af286bca1b), UWORD(0xbca1af286bca1af2), UWORD(0xe50d79435e50d794), UWORD(0xd79435e50d79435) },  // 19
    { UWORD(0xd37a6f4de9bd37a7), UWORD(0x4de9bd37a6f4de9b), UWORD(0x42c8590b21642c85), UWORD(0xb21642c8590b216) },  // 23
    { UWORD(0x34f72c234f72c235), UWORD(0xc234f72c234f72c2), UWORD(0xd3dcb08d3dcb08d3), UWORD(0x8d3dcb08d3dcb08) },  // 29
    { UWORD(0xef7bdef7bdef7bdf), UWORD(0xdef7bdef7bdef7bd), UWORD(0x8421084210842108), UWORD(0x842108421084210) },  // 31
    { UWORD(0x14c1bacf914c1bad), UWORD(0xc1bacf914c1bacf9), UWORD(0x5306eb3e45306eb3), UWORD(0x6eb3e45306eb3e4) },  // 37
    { UWORD(0x8f9c18f9c18f9c19), UWORD(0x18f9c18f9c18f9c1), UWORD(0x63e7063e7063e706), UWORD(0x63e7063e7063e70) },  // 41
    { UWORD(0x82fa0be82fa0be83), UWORD(0xbe82fa0be82fa0be), UWORD(0xf417d05f417d05f4), UWORD(0x5f417d05f417d05) },  // 43
    { UWORD(0x51b3bea3677d46cf), UWORD(0x3677d46cefa8d9df), UWORD(0x882b9310572620ae), UWORD(0x572620ae4c415c9) },  // 47
    { UWORD(0x21cfb2b78c13521d), UWORD(0x13521cfb2b78c135), UWORD(0x4873ecade304d487), UWORD(0x4d4873ecade304d) },  // 53
    { UWORD(0xcbeea4e1a08ad8f3), UWORD(0x8f2fba9386822b63), UWORD(0x15b1e5f75270d045), UWORD(0x456c797dd49c341) },  // 59
    { UWORD(0x4fbcda3ac10c9715), UWORD(0x14fbcda3ac10c971), UWORD(0x4325c53ef368eb04), UWORD(0x4325c53ef368eb0) },  // 61
    { UWORD(0xf0b7672a07a44c6b), UWORD(0xc2dd9ca81e9131ab), UWORD(0x40f4898d5f85bb39), UWORD(0x3d226357e16ece5) },  // 67
    { UWORD(0x193d4bb7e327a977), UWORD(0x4f52edf8c9ea5dbf), UWORD(0x240e6c2b4481cd85), UWORD(0x39b0ad12073615a) },  // 71
    { UWORD(0x7e3f1f8fc7e3f1f9), UWORD(0x3f1f8fc7e3f1f8fc), UWORD(0x70381c0e070381c), UWORD(0x381c0e070381c0e) },  // 73
    { UWORD(0x9b8b577e613716af), UWORD(0xd5df984dc5abbf30), UWORD(0xa5440cf6474a8819), UWORD(0x33d91d2a2067b23) },  // 79
    { UWORD(0xa3784a062b2e43db), UWORD(0x2818acb90f6bf3a9), UWORD(0x6f0940c565c87b5f), UWORD(0x3159721ed7e7534) },  // 83
    { UWORD(0xf47e8fd1fa3f47e9), UWORD(0xd1fa3f47e8fd1fa3), UWORD(0xc0b81702e05c0b81), UWORD(0x2e05c0b81702e05) },  // 89
    { UWORD(0xa3a0fd5c5f02a3a1), UWORD(0x5f02a3a0fd5c5f02), UWORD(0xa0fd5c5f02a3a0fd), UWORD(0x2a3a0fd5c5f02a3) },  // 97
    { UWORD(0x3a4c0a237c32b16d), UWORD(0xc32b16cfd7720f35), UWORD(0xc83cd4e930288df0), UWORD(0x288df0cac5b3f5d) },  // 101
    { UWORD(0xdab7ec1dd3431b57), UWORD(0xd0c6d5bf60ee9a18), UWORD(0x88b2f392a409f116), UWORD(0x27c45979c95204f) },  // 103
    { UWORD(0x77a04c8f8d28ac43), UWORD(0xa2b10bf66e0e5aea), UWORD(0xdc1cb5d4ef40991f), UWORD(0x2647c69456217ec) },  // 107
    { UWORD(0xa6c0964fda6c0965), UWORD(0xc0964fda6c0964fd), UWORD(0x9b02593f69b02593), UWORD(0x2593f69b02593f6) },  // 109
    { UWORD(0x90fdbc090fdbc091), UWORD(0xc090fdbc090fdbc0), UWORD(0x43f6f0243f6f0243), UWORD(0x243f6f0243f6f02) },  // 113
    { UWORD(0x7efdfbf7efdfbf7f), UWORD(0xbf7efdfbf7efdfbf), UWORD(0x408102040810204), UWORD(0x204081020408102) },  // 127
    { UWORD(0x3e88cb3c9484e2b), UWORD(0xf82ee6986d6f63aa), UWORD(0x7f05dcd30dadec75), UWORD(0x1f44659e4a42715) },  // 131
    { UWORD(0xe21a291c077975b9), UWORD(0x21a291c077975b8f), UWORD(0x701de5d6e3f8868a), UWORD(0x1de5d6e3f8868a4) },  // 137
    { UWORD(0x3aef6ca970586723), UWORD(0xa2126ad1f4f31ba0), UWORD(0x17f14424d5a3e9e6), UWORD(0x1d77b654b82c339) },  // 139
    { UWORD(0xdf5b0f768ce2cabd), UWORD(0x93c225cc74d50c06), UWORD(0xaf3f920a4f089731), UWORD(0x1b7d6c3dda338b2) },  // 149
    { UWORD(0x6fe4dfc9bf937f27), UWORD(0x26fe4dfc9bf937f2), UWORD(0x1b2036406c80d901), UWORD(0x1b2036406c80d90) },  // 151
    { UWORD(0x5b4fe5e92c0685b5), UWORD(0x685b4fe5e92c068), UWORD(0x16d3f97a4b01a16d), UWORD(0x1a16d3f97a4b01a) },  // 157
    { UWORD(0x1f693a1c451ab30b), UWORD(0x8bc775ca99ea0324), UWORD(0x59857f36f825b178), UWORD(0x1920fb49d0e228d) },  // 163
    { UWORD(0x8d07aa27db35a717), UWORD(0x513ed9ad38b7f3bc), UWORD(0x4b1d20310dcbe157), UWORD(0x1886e5f0abb0499) },  // 167
    { UWORD(0x882383b30d516325), UWORD(0x133caba736c05eb4), UWORD(0x458c93fa14b77dc7), UWORD(0x17ad2208e0ecc35) },  // 173
    { UWORD(0xed6866f8d962ae7b), UWORD(0xe4d3aa30a02dc3e), UWORD(0xb1573d7f48f044a5), UWORD(0x16e1f76b4337c6c) },  // 179
    { UWORD(0x3454dca410f8ed9d), UWORD(0x6fbc1c498c05a84f), UWORD(0x3e3b673fa57b0cba), UWORD(0x16a13cd15372904) },  // 181
    { UWORD(0x1d7ca632ee936f3f), UWORD(0x7749b79f7f547096), UWORD(0x22d9218202ae3da7), UWORD(0x1571ed3c506b39a) },  // 191
    { UWORD(0x70bf015390948f41), UWORD(0x90948f40feac6f6b), UWORD(0x6f6b70bf01539094), UWORD(0x15390948f40feac) },  // 193
    { UWORD(0xc96bdb9d3d137e0d), UWORD(0xbb207cc0532ae21), UWORD(0x4f44df833facd51d), UWORD(0x14cab88725af6e7) },  // 197
    { UWORD(0x2697cc8aef46c0f7), UWORD(0x7a3607b7f5b5630e), UWORD(0xa21727e120292a73), UWORD(0x149539e3b2d066e) },  // 199
    { UWORD(0xc0e8f2a76e68575b), UWORD(0x2f514a026d31be7b), UWORD(0x53b7342bad7f64b3), UWORD(0x13698df3de07479) },  // 211
    { UWORD(0x687763dfdb43bb1f), UWORD(0xdd8f7f6d0eec7bfb), UWORD(0x3840497889c2024b), UWORD(0x125e22708092f11) },  // 223
    { UWORD(0x1b10ea929ba144cb), UWORD(0x766a024168e18cf8), UWORD(0x75494dd0a2657f6f), UWORD(0x120b470c67c0d88) },  // 227
    { UWORD(0x1d10c4c0478bbced), UWORD(0xc4c0478bbcecfee), UWORD(0x313011e2ef3b3fb8), UWORD(0x11e2ef3b3fb8744) },  // 229
    { UWORD(0x63fb9aeb1fdcd759), UWORD(0x758fee6bac7f735d), UWORD(0x46514e02328a7011), UWORD(0x119453808ca29c0) },  // 233
    { UWORD(0x64afaa4f437b2e0f), UWORD(0x77f76e538c5167e), UWORD(0xa0ab617909a3e202), UWORD(0x112358e75d30336) },  // 239
    { UWORD(0xf010fef010fef011), UWORD(0x10fef010fef010fe), UWORD(0xef010fef010fef01), UWORD(0x10fef010fef010f) },  // 241
    { UWORD(0x28cbfbeb9a020a33), UWORD(0xa020a32fefae6808), UWORD(0x465fdf5cd0105197), UWORD(0x105197f7d734041) },  // 251
    { UWORD(0xff00ff00ff00ff01), UWORD(0xff00ff00ff00ff00), UWORD(0xff00ff00ff00ff), UWORD(0xff00ff00ff00ff) },  // 257
    { UWORD(0xd624fd1470e99cb7), UWORD(0xf836826ef73d52bc), UWORD(0x653b605d71e2cc69), UWORD(0xf92fb2211855a8) },  // 263
    { UWORD(0x8fb3ddbd6205b5c5), UWORD(0x3ce8354b2ea1c8cd), UWORD(0x363ecf76f58816d7), UWORD(0xf3a0d52cba8723) },  // 269
    { UWORD(0xd57da36ca27acdef), UWORD(0x8715ba188f963302), UWORD(0xfa5504b926bb0a64), UWORD(0xf1d48bcee0d399) },  // 271
    { UWORD(0xee70c03b25e4463d), UWORD(0xb25e4463cff13686), UWORD(0xa1bb9c300ec97911), UWORD(0xec979118f3fc4d) },  // 277
    { UWORD(0xc5b1a6b80749cb29), UWORD(0x6c69ae01d272ca3f), UWORD(0x5c03a4e5947f8b63), UWORD(0xe939651fe2d8d3) },  // 281
    { UWORD(0x47768073c9b97113), UWORD(0xf26e5c44bfc61b23), UWORD(0xd91a3bb4039e4dcb), UWORD(0xe79372e225fe30) },  // 283
    { UWORD(0x2591e94884ce32ad), UWORD(0xb07dd0d1b15d7cf1), UWORD(0x5f3c49647a522133), UWORD(0xdfac1f74346c57) },  // 293
    { UWORD(0xf02806abc74be1fb), UWORD(0xd2f87ebfcaa1c5a0), UWORD(0x50e2d078140355e3), UWORD(0xd578e97c3f5fe5) },  // 307
    { UWORD(0x7ec3e8f3a7198487), UWORD(0xbe25dd6d7aa646ca), UWORD(0xab3726b02782e18b), UWORD(0xd2ba083b445250) },  // 311
    { UWORD(0x58550f8a39409d09), UWORD(0xbc1d71afd8bdc034), UWORD(0x7423fcba7aaf075c), UWORD(0xd161543e28e502) },  // 313
};

/* Let d be the largest trial prime. If x < 2^(2*FLINT_BITS) - d*2^FLINT_BITS, then
   a < comparison of the high limbs in the Granlund-Montgomery divisibility test
   is guaranteed to give the same result as a <= comparison of both limbs. */
#define FAST_TRIAL_BOUND (UWORD_MAX - 313)

/* The fallback full comparison will rarely be used except e.g. when
   generating decreasing primes from 2^128. */
static int
n_ll_le(ulong ahi, ulong alo, ulong bhi, ulong blo)
{
#ifdef __SIZEOF_INT128__
    return (((__uint128_t) ahi << 64) | alo) <=
           (((__uint128_t) bhi << 64) | blo);
#else
    return (ahi < bhi) | ((ahi == bhi) & (alo <= blo));
#endif

}


/* FIXME: until we have good double-limb implementation, mock one up. */
#include <fmpz.h>

static int
n_ll_fallback_is_prime_base2(ulong nhi, ulong nlo)
{
    int res;
    fmpz_t nn;
    fmpz_init(nn);
    fmpz_set_uiui(nn, nhi, nlo);
    fmpz base2 = 2;
    if (!fmpz_is_strong_probabprime(nn, &base2))
        res = 0;
    else
        res = -1;
    fmpz_clear(nn);
    return res;
}

int
n_ll_is_prime(ulong nhi, ulong nlo)
{
    ulong n1 = nhi, n0 = nlo;
    ulong a0, a1;
    ulong b0, b1;
    slong i;

    FLINT_ASSERT(nhi != 0);

    /* Trial division by 2 */
    if (n0 % 2 == 0)
        return 0;

    /* Granlund-Montgomery trial division by 3, 5, ..., 137.
       On Zen 3, unrolling by two and combining comparisons bitwise
       is faster than a more branchy version. */
#define ainv0 ll_trial_primes[i][0]
#define ainv1 ll_trial_primes[i][1]
#define abound0 ll_trial_primes[i][2]
#define abound1 ll_trial_primes[i][3]
#define binv0 ll_trial_primes[i + 1][0]
#define binv1 ll_trial_primes[i + 1][1]
#define bbound0 ll_trial_primes[i + 1][2]
#define bbound1 ll_trial_primes[i + 1][3]

    if (nhi < FAST_TRIAL_BOUND)
    {
        for (i = 0; i < 32; i += 2)
        {
            mullo_2x2(&a1, &a0, n1, n0, ainv1, ainv0);
            mullo_2x2(&b1, &b0, n1, n0, binv1, binv0);
            if ((a1 < abound1) | (b1 < bbound1))
                return 0;
        }
    }
    else
    {
        for (i = 0; i < 32; i += 2)
        {
            mullo_2x2(&a1, &a0, n1, n0, ainv1, ainv0);
            mullo_2x2(&b1, &b0, n1, n0, binv1, binv0);
            if (n_ll_le(a1, a0, abound1, abound0) | n_ll_le(b1, b0, bbound1, bbound0))
                return 0;
        }
    }

    if (n1 > SWbound[1] || (n1 == SWbound[1] && n0 >= SWbound[0]))
    {
        /* 6 * 2 * (357913941 * 2^64) < 2^192, so we can use the fast
           base-2 Miller-Rabin code */
        if (n1 < UWORD(357913941))
            return n_ll_small_is_prime_mr(n1, n0, 1);
        else /* The base-2 test is slower now, so do some extra trial division. */
        {
            if (nhi < FAST_TRIAL_BOUND)
            {
                for (i = 32; i < 64; i += 2)
                {
                    mullo_2x2(&a1, &a0, n1, n0, ainv1, ainv0);
                    mullo_2x2(&b1, &b0, n1, n0, binv1, binv0);
                    if ((a1 < abound1) | (b1 < bbound1))
                        return 0;
                }
            }
            else
            {
                for (i = 32; i < 64; i += 2)
                {
                    mullo_2x2(&a1, &a0, n1, n0, ainv1, ainv0);
                    mullo_2x2(&b1, &b0, n1, n0, binv1, binv0);
                    if (n_ll_le(a1, a0, abound1, abound0) | n_ll_le(b1, b0, bbound1, bbound0))
                        return 0;
                }
            }

            return n_ll_fallback_is_prime_base2(n1, n0);
        }
    }

    return n_ll_small_is_prime_mr(n1, n0, 0);
}

#endif

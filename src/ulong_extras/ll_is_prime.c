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

#if FLINT_BITS == 32

/* Todo: could implement a reasonable 64-bit primality test here. */
int
n_ll_is_prime(ulong nhi, ulong nlo)
{
    return -1;
}

#else

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
    ulong x1, ulong x0, ulong m1, ulong m0, ulong minv1, ulong minv0)
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

    /* Do the first two rounds without modular reduction, since the result
       certainly fits in one word. To do: could do one more round, with a
       two-word result. */
    if (ebits >= 3)
    {
        ulong bit2, bit1;

        bit1 = exp[(ebits - 2) / 64] >> ((ebits - 2) % 64);
        bit2 = exp[(ebits - 3) / 64] >> ((ebits - 3) % 64);

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
        elimbs = (ebits + 64) / 64;
    }

    if (ebits >= 64)
        ebits -= 64;

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

    if (ebits >= 6)
    {
        ulong exp_top;
        ulong shift = ebits - 6;

        if (shift >= 64)
            exp_top = exp[1] >> (shift - 64);
        else if (ebits >= 64)
            exp_top = (exp[1] << (64 - shift)) | (exp[0] >> shift);
        else
            exp_top = exp[0] >> shift;

        xr0 = UWORD(1) << exp_top;
        ebits -= 6;
        elimbs = (ebits + 64) / 64;
        if (ebits >= 64)
            ebits -= 64;
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
n_ll_small_is_prime_mr_13base(ulong n1, ulong n0)
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

/* Low and high word of 1/p mod 2^128; low and high word of 2^128 / p */
static const ulong ll_trial_primes[302][4] = {
    { UWORD(0xaaaaaaaaaaaaaaab), UWORD(0xaaaaaaaaaaaaaaaa), UWORD(0xaaaaaaaaaaaaaaaa), UWORD(0x2aaaaaaaaaaaaaaa) },  // 3
    { UWORD(0xcccccccccccccccd), UWORD(0xcccccccccccccccc), UWORD(0x9999999999999999), UWORD(0x1999999999999999) },  // 5
    { UWORD(0x6db6db6db6db6db7), UWORD(0xb6db6db6db6db6db), UWORD(0x2492492492492492), UWORD(0x1249249249249249) },  // 7
    { UWORD(0x2e8ba2e8ba2e8ba3), UWORD(0xa2e8ba2e8ba2e8ba), UWORD(0xba2e8ba2e8ba2e8b), UWORD(0xba2e8ba2e8ba2e8) },  // 11
    { UWORD(0x4ec4ec4ec4ec4ec5), UWORD(0xc4ec4ec4ec4ec4ec), UWORD(0x9d89d89d89d89d89), UWORD(0x9d89d89d89d89d8) },  // 13
    { UWORD(0xf0f0f0f0f0f0f0f1), UWORD(0xf0f0f0f0f0f0f0f0), UWORD(0x8787878787878787), UWORD(0x787878787878787) },  // 17
    { UWORD(0x86bca1af286bca1b), UWORD(0xbca1af286bca1af2), UWORD(0xf286bca1af286bca), UWORD(0x6bca1af286bca1a) },  // 19
    { UWORD(0xd37a6f4de9bd37a7), UWORD(0x4de9bd37a6f4de9b), UWORD(0x21642c8590b21642), UWORD(0x590b21642c8590b) },  // 23
    { UWORD(0x34f72c234f72c235), UWORD(0xc234f72c234f72c2), UWORD(0x69ee58469ee58469), UWORD(0x469ee58469ee584) },  // 29
    { UWORD(0xef7bdef7bdef7bdf), UWORD(0xdef7bdef7bdef7bd), UWORD(0x4210842108421084), UWORD(0x421084210842108) },  // 31
    { UWORD(0x14c1bacf914c1bad), UWORD(0xc1bacf914c1bacf9), UWORD(0x2983759f22983759), UWORD(0x3759f22983759f2) },  // 37
    { UWORD(0x8f9c18f9c18f9c19), UWORD(0x18f9c18f9c18f9c1), UWORD(0x31f3831f3831f383), UWORD(0x31f3831f3831f38) },  // 41
    { UWORD(0x82fa0be82fa0be83), UWORD(0xbe82fa0be82fa0be), UWORD(0xfa0be82fa0be82fa), UWORD(0x2fa0be82fa0be82) },  // 43
    { UWORD(0x51b3bea3677d46cf), UWORD(0x3677d46cefa8d9df), UWORD(0xc415c9882b931057), UWORD(0x2b9310572620ae4) },  // 47
    { UWORD(0x21cfb2b78c13521d), UWORD(0x13521cfb2b78c135), UWORD(0xa439f656f1826a43), UWORD(0x26a439f656f1826) },  // 53
    { UWORD(0xcbeea4e1a08ad8f3), UWORD(0x8f2fba9386822b63), UWORD(0x8ad8f2fba9386822), UWORD(0x22b63cbeea4e1a0) },  // 59
    { UWORD(0x4fbcda3ac10c9715), UWORD(0x14fbcda3ac10c971), UWORD(0x2192e29f79b47582), UWORD(0x2192e29f79b4758) },  // 61
    { UWORD(0xf0b7672a07a44c6b), UWORD(0xc2dd9ca81e9131ab), UWORD(0xa07a44c6afc2dd9c), UWORD(0x1e9131abf0b7672) },  // 67
    { UWORD(0x193d4bb7e327a977), UWORD(0x4f52edf8c9ea5dbf), UWORD(0x12073615a240e6c2), UWORD(0x1cd85689039b0ad) },  // 71
    { UWORD(0x7e3f1f8fc7e3f1f9), UWORD(0x3f1f8fc7e3f1f8fc), UWORD(0x381c0e070381c0e), UWORD(0x1c0e070381c0e07) },  // 73
    { UWORD(0x9b8b577e613716af), UWORD(0xd5df984dc5abbf30), UWORD(0xd2a2067b23a5440c), UWORD(0x19ec8e951033d91) },  // 79
    { UWORD(0xa3784a062b2e43db), UWORD(0x2818acb90f6bf3a9), UWORD(0x3784a062b2e43daf), UWORD(0x18acb90f6bf3a9a) },  // 83
    { UWORD(0xf47e8fd1fa3f47e9), UWORD(0xd1fa3f47e8fd1fa3), UWORD(0xe05c0b81702e05c0), UWORD(0x1702e05c0b81702) },  // 89
    { UWORD(0xa3a0fd5c5f02a3a1), UWORD(0x5f02a3a0fd5c5f02), UWORD(0xd07eae2f8151d07e), UWORD(0x151d07eae2f8151) },  // 97
    { UWORD(0x3a4c0a237c32b16d), UWORD(0xc32b16cfd7720f35), UWORD(0xe41e6a74981446f8), UWORD(0x1446f86562d9fae) },  // 101
    { UWORD(0xdab7ec1dd3431b57), UWORD(0xd0c6d5bf60ee9a18), UWORD(0xc45979c95204f88b), UWORD(0x13e22cbce4a9027) },  // 103
    { UWORD(0x77a04c8f8d28ac43), UWORD(0xa2b10bf66e0e5aea), UWORD(0x6e0e5aea77a04c8f), UWORD(0x1323e34a2b10bf6) },  // 107
    { UWORD(0xa6c0964fda6c0965), UWORD(0xc0964fda6c0964fd), UWORD(0x4d812c9fb4d812c9), UWORD(0x12c9fb4d812c9fb) },  // 109
    { UWORD(0x90fdbc090fdbc091), UWORD(0xc090fdbc090fdbc0), UWORD(0x21fb78121fb78121), UWORD(0x121fb78121fb781) },  // 113
    { UWORD(0x7efdfbf7efdfbf7f), UWORD(0xbf7efdfbf7efdfbf), UWORD(0x204081020408102), UWORD(0x102040810204081) },  // 127
    { UWORD(0x3e88cb3c9484e2b), UWORD(0xf82ee6986d6f63aa), UWORD(0xbf82ee6986d6f63a), UWORD(0xfa232cf252138a) },  // 131
    { UWORD(0xe21a291c077975b9), UWORD(0x21a291c077975b8f), UWORD(0x380ef2eb71fc4345), UWORD(0xef2eb71fc43452) },  // 137
};

int
n_ll_is_prime(ulong nhi, ulong nlo)
{
    ulong n1 = nhi, n0 = nlo;
    ulong t1, t0;
    slong i;

    FLINT_ASSERT(nhi != 0);

    /* Trial division by 2 */
    if (n0 % 2 == 0)
        return 0;

    /* To do: move this check below trial division. Currently this
       abort is done early to avoid duplicate trial division in
       fmpz_is_prime. */
    if (n1 > SWbound[1] || (n1 == SWbound[1] && n0 >= SWbound[0]))
        return -1;

    /* Granlund-Montgomery trial division by 3, 5, ..., 137.
       On Zen 3, unrolling by two and combining comparisons bitwise
       is faster than a more branchy version. */
    for (i = 0; i < 32; i += 2)
    {
        ulong v1, v0, u1, u0;
        ulong w1, w0, q1, q0, r1, r0;

        u0 = ll_trial_primes[i][0];
        u1 = ll_trial_primes[i][1];
        v0 = ll_trial_primes[i][2];
        v1 = ll_trial_primes[i][3];

        q0 = ll_trial_primes[i + 1][0];
        q1 = ll_trial_primes[i + 1][1];
        w0 = ll_trial_primes[i + 1][2];
        w1 = ll_trial_primes[i + 1][3];

        mullo_2x2(&t1, &t0, n1, n0, u1, u0);
        mullo_2x2(&r1, &r0, n1, n0, q1, q0);

        if ((t1 < v1) | ((t1 == v1) & (t0 <= v0)) | (r1 < v1) | ((r1 == w1) & (r0 <= w0)))
            return 0;
    }

    return n_ll_small_is_prime_mr_13base(n1, n0);
}

#endif

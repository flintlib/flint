/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPN_EXTRAS_RECIP_SMALL_H
#define MPN_EXTRAS_RECIP_SMALL_H

/* Reciprocals at four to eleven limbs by Newton-type steps on short
   products, as inline functions meant for a compile-time size: the
   engine of the small-size approximate divisions (divapprox_small.c)
   and of the dense-table kernels of mp_real.  64-bit limbs only. */

#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"

#if FLINT_BITS == 64

#define FLINT_MPN_RECIP_SMALL_MAX 11

/* r = u +- v over a compile-time number of limbs 2 .. 11, no carry out */
#define _FLINT_MPN_RECIP_NN_ADD_9(r, u, v) add_sssssssssaaaaaaaaaaaaaaaaaa((r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN_SUB_9(r, u, v) sub_dddddddddmmmmmmmmmsssssssss((r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN_ADD_10(r, u, v) add_ssssssssssaaaaaaaaaaaaaaaaaaaa((r)[9], (r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[9], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[9], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN_SUB_10(r, u, v) sub_ddddddddddmmmmmmmmmmssssssssss((r)[9], (r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[9], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[9], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN_ADD_11(r, u, v) add_sssssssssssaaaaaaaaaaaaaaaaaaaaaa((r)[10], (r)[9], (r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[10], (u)[9], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[10], (v)[9], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN_SUB_11(r, u, v) sub_dddddddddddmmmmmmmmmmmsssssssssss((r)[10], (r)[9], (r)[8], (r)[7], (r)[6], (r)[5], (r)[4], (r)[3], (r)[2], (r)[1], (r)[0], (u)[10], (u)[9], (u)[8], (u)[7], (u)[6], (u)[5], (u)[4], (u)[3], (u)[2], (u)[1], (u)[0], (v)[10], (v)[9], (v)[8], (v)[7], (v)[6], (v)[5], (v)[4], (v)[3], (v)[2], (v)[1], (v)[0])
#define _FLINT_MPN_RECIP_NN(OP, n, r, u, v) \
    do { switch (n) { \
        case 2: NN_##OP##_2(r, u, v); break; \
        case 3: NN_##OP##_3(r, u, v); break; \
        case 4: NN_##OP##_4(r, u, v); break; \
        case 5: NN_##OP##_5(r, u, v); break; \
        case 6: NN_##OP##_6(r, u, v); break; \
        case 7: NN_##OP##_7(r, u, v); break; \
        case 8: NN_##OP##_8(r, u, v); break; \
        case 9: _FLINT_MPN_RECIP_NN_##OP##_9(r, u, v); break; \
        case 10: _FLINT_MPN_RECIP_NN_##OP##_10(r, u, v); break; \
        case 11: _FLINT_MPN_RECIP_NN_##OP##_11(r, u, v); break; \
        default: FLINT_ASSERT(0); FLINT_UNREACHABLE; } } while (0)

/* One Newton-type step at n limbs (4 <= n <= FLINT_MPN_RECIP_SMALL_MAX)
   from (Iq, q) <= 1/X, X = 1 + F: with E = 1 - X Iq (formed from the
   exact product of F and Iq, its low q limbs dropped: at most one ulp
   high), whose top z = 2 limbs vanish (Iq is within 2^-128 of 1/X),
   I = (Iq - lower ulps) + Iq (E + E^2 + ... + E^(order - 1)), the
   products of E' = E B^z (n - z limbs) truncated (under one ulp each).
   With d = 1/X - Iq, the true error of the step is below
   X^(order - 1) d^order, so that I lies within that, plus lower + 3
   ulps, below 1/X, and at least lower - 1 ulps below it (up to the
   true error). */
FLINT_FORCE_INLINE void
_flint_mpn_recip_small_step(nn_ptr I, nn_srcptr F, slong n, nn_srcptr Iq,
    slong q, int order, ulong lower)
{
    const slong z = 2;
    ulong P[2 * FLINT_MPN_RECIP_SMALL_MAX + 2], E[FLINT_MPN_RECIP_SMALL_MAX + 1];
    ulong Z[FLINT_MPN_RECIP_SMALL_MAX + 1], S[FLINT_MPN_RECIP_SMALL_MAX + 1];
    slong i, m;

    /* the contract, also as hints for the compilers' range analysis
       (which otherwise warns about indices of sizes never used) */
    FLINT_ASSERT(n >= 4 && n <= FLINT_MPN_RECIP_SMALL_MAX);
    FLINT_ASSERT(q >= 1 && q <= n);
    if (n < 4 || n > FLINT_MPN_RECIP_SMALL_MAX || q < 1 || q > n)
        FLINT_UNREACHABLE;

    flint_mpn_mul(P, F, n, Iq, q);
    for (i = 0; i < n; i++)
        Z[i] = 0;
    for (i = 0; i < n - q; i++)
        E[i] = 0;
    for (i = 0; i < q; i++)
        E[n - q + i] = Iq[i];
    _FLINT_MPN_RECIP_NN(ADD, n, E, E, P + q);
    _FLINT_MPN_RECIP_NN(SUB, n, E, Z, E);
#if FLINT_WANT_ASSERT
    for (i = n - z; i < n; i++)
        FLINT_ASSERT(E[i] == 0);
#endif

    if (order >= 3)
    {
        /* E^2 = E'^2 B^-2z: limb k of the high product S of E' (n - z
           limbs) goes to the position k - z */
        flint_mpn_sqrhigh_n(S, E, n - z);
        if (order >= 4)
        {
            /* E^3 = E'^3 B^-3z: limb k of the high product of the top
               m = n - 2z limbs of S and E' goes to the position k - z */
            m = n - 2 * z;
            FLINT_ASSERT(m >= 1);
            if (m < 1)
                FLINT_UNREACHABLE;
            flint_mpn_mulhigh_n(P, S + n - z - m, E + n - z - m, m);
            for (i = m; i < n; i++)
                P[i] = 0;
            _FLINT_MPN_RECIP_NN(ADD, n - z, E, E, P + z);
        }
        for (i = n - z; i < n; i++)
            S[i] = 0;
        _FLINT_MPN_RECIP_NN(ADD, n - z, E, E, S + z);
    }

    /* the product of E' and Iq at weight B^-(n + q), its limbs
       q .. n - z + q - 1 at the positions 0 .. n - z - 1 */
    if (n - z >= q)
        flint_mpn_mul(P, E, n - z, Iq, q);
    else
        flint_mpn_mul(P, Iq, q, E, n - z);
    for (i = n - z + q; i < n + q; i++)
        P[i] = 0;
    for (i = 0; i < n - q; i++)
        I[i] = 0;
    for (i = 0; i < q; i++)
        I[n - q + i] = Iq[i];
    Z[0] = lower;
    _FLINT_MPN_RECIP_NN(SUB, n, I, I, Z);
    _FLINT_MPN_RECIP_NN(ADD, n, I, I, P + q);
}

/* (I, n) <= 1/(1 + F) for the n-limb fraction F (4 <= n <= FLINT_MPN_RECIP_SMALL_MAX), within
   six ulps of B^-n -- or within 2^-bits when the caller needs fewer
   bits, which may lower the order of the last step.

   The one-limb reciprocal R of the top bits (I_0 = R/B, within 2^-61.9
   below 1/X, X = 1 + F); then a cubic step at three limbs on products
   by the single limb R, I_1 = I_0 (1 + E + E^2), E = 1 - X I_0 (the
   true error I_0 E^3 / (1 - E) below 2^-183; E is formed at most two
   ulps high, so I_1 is lowered by three ulps of B^-3: within 2^-182.9
   below 1/X); then one step at n limbs on products by the three limbs
   of I_1 (E < 2^-181.9 has two zero top limbs), of order 2 (true error
   below 2^-364), 3 (2^-546) or 4 (2^-728) as n needs, lowered by two
   ulps (E is formed at most one ulp high; the up to three truncated
   products lose under one ulp each).  So I lies in (1/X - 5 - 2^-100,
   1/X - 1 + 2^-100] ulps. */
FLINT_FORCE_INLINE void
_flint_mpn_recip_small(nn_ptr I, nn_srcptr F, slong n, slong bits)
{
    ulong R, rr, I1[3], E1[3], G1[3];
    ulong p3, p2, p1, p0, q3, q2, q1;

    FLINT_ASSERT(n >= 4 && n <= FLINT_MPN_RECIP_SMALL_MAX);

    /* R / B <= 1/X: with Xt = B/2 + floor(F_top / 2) <= X B/2,
       (B^2/2 - 1) / Xt - 2 < B / X, and the floor loses at most three
       units more */
    udiv_qrnnd(R, rr, (UWORD(1) << (FLINT_BITS - 1)) - 1, ~UWORD(0),
        (UWORD(1) << (FLINT_BITS - 1)) | (F[n - 1] >> 1));
    (void) rr;
    R -= 2;

    /* three limbs: E = 1 - (I_0 + F_3 I_0) mod B^3 (the low word of
       F_3 I_0 dropped; I_0 + F_3 I_0 <= X I_0 <= 1, equality wrapping
       to E = 0) */
    umul_ppmm(p1, p0, F[n - 3], R);
    (void) p0;
    umul_ppmm(q2, q1, F[n - 2], R);
    add_ssaaaa(p2, p1, UWORD(0), p1, q2, q1);
    umul_ppmm(q3, q2, F[n - 1], R);
    add_ssaaaa(p3, p2, UWORD(0), p2, q3, q2);
    p3 += R;
    sub_dddmmmsss(E1[2], E1[1], E1[0], UWORD(0), UWORD(0), UWORD(0), p3, p2, p1);

    /* G = E + E^2 */
    flint_mpn_sqrhigh_n(G1, E1, 3);
    NN_ADD_3(G1, G1, E1);

    /* I_1 = (I_0 - 3 ulps) + I_0 G (the low word dropped), below 1 */
    umul_ppmm(p1, p0, G1[0], R);
    (void) p0;
    umul_ppmm(q2, q1, G1[1], R);
    add_ssaaaa(p2, p1, UWORD(0), p1, q2, q1);
    umul_ppmm(q3, q2, G1[2], R);
    add_ssaaaa(p3, p2, UWORD(0), p2, q3, q2);
    add_sssaaaaaa(I1[2], I1[1], I1[0], R - 1, ~UWORD(0), ~UWORD(0) - 2, p3, p2, p1);

    /* one step at n limbs: quadratic, cubic from 6 limbs (E^2 >
       2^-364), quartic from 8 (E^3 > 2^-546) */
    _flint_mpn_recip_small_step(I, F, n, I1, 3,
        (n >= 8 && bits > 540) ? 4 : (n >= 6 && bits > 360) ? 3 : 2, 2);
}

#endif

#endif

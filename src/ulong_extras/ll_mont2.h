/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Opus 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
   Two-word Montgomery arithmetic, shared by the double-word factoring
   routines (ll_factor_rho.c, ll_factor_ecm.c).

   The modulus n = {n1,n0} must be odd with n1 != 0, and np = -n^(-1) mod
   2^FLINT_BITS.  Values are held as ordinary two-word residues in [0,n);
   "Montgomery form" of a is a*R mod n with R = 2^(2*FLINT_BITS).

   Note that gcd(a*R mod n, n) = gcd(a, n) because R is a unit mod n, so
   factor detection needs no conversion out of Montgomery form.
*/

#ifndef LL_MONT2_H
#define LL_MONT2_H

#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

#if FLINT_BITS == 64

/*
   Products are computed with the existing mpn_extras macros:
   FLINT_MPN_MUL_2X2 for a general 2x2 product and FLINT_MPN_SQR_2X2 for a
   square, which saves one multiplication.
*/

/*
   Montgomery reduction of {p3,p2,p1,p0} modulo {n1,n0}.  A fifth
   accumulator word is carried because the intermediate sum can exceed
   four words.  The result is reduced to [0, n).
*/
#define LL_REDC_2(r1, r0, p3, p2, p1, p0, n1, n0, np)                       \
    do {                                                                    \
        ulong __T4 = 0, __T3 = (p3), __T2 = (p2), __T1 = (p1), __T0 = (p0); \
        ulong __m, __h, __l, __cy;                                          \
                                                                            \
        __m = __T0 * (np);                                                  \
        umul_ppmm(__h, __l, __m, (n0));                                     \
        add_ssaaaa(__cy, __T0, 0, __T0, __h, __l);   /* __T0 becomes 0 */   \
        umul_ppmm(__h, __l, __m, (n1));                                     \
        add_ssaaaa(__h, __l, __h, __l, 0, __cy);                            \
        add_ssssaaaaaaaa(__T4, __T3, __T2, __T1,                            \
                         __T4, __T3, __T2, __T1, 0, 0, __h, __l);           \
                                                                            \
        __m = __T1 * (np);                                                  \
        umul_ppmm(__h, __l, __m, (n0));                                     \
        add_ssaaaa(__cy, __T1, 0, __T1, __h, __l);   /* __T1 becomes 0 */   \
        umul_ppmm(__h, __l, __m, (n1));                                     \
        add_ssaaaa(__h, __l, __h, __l, 0, __cy);                            \
        add_sssaaaaaa(__T4, __T3, __T2, __T4, __T3, __T2, 0, __h, __l);     \
                                                                            \
        if (__T4 || __T3 > (n1) || (__T3 == (n1) && __T2 >= (n0)))          \
            sub_ddmmss(__T3, __T2, __T3, __T2, (n1), (n0));                 \
                                                                            \
        (r1) = __T3;                                                        \
        (r0) = __T2;                                                        \
    } while (0)

#define LL_MULMOD(r1, r0, a1, a0, b1, b0, n1, n0, np)                       \
    do {                                                                    \
        ulong __q3, __q2, __q1, __q0;                                       \
        FLINT_MPN_MUL_2X2(__q3, __q2, __q1, __q0, a1, a0, b1, b0);          \
        LL_REDC_2(r1, r0, __q3, __q2, __q1, __q0, n1, n0, np);              \
    } while (0)

/* as LL_MULMOD, but squaring {a1,a0} */
#define LL_SQRMOD(r1, r0, a1, a0, n1, n0, np)                               \
    do {                                                                    \
        ulong __q3, __q2, __q1, __q0;                                       \
        FLINT_MPN_SQR_2X2(__q3, __q2, __q1, __q0, a1, a0);                  \
        LL_REDC_2(r1, r0, __q3, __q2, __q1, __q0, n1, n0, np);              \
    } while (0)

/* add {b1,b0} to {a1,a0} modulo {n1,n0}; inputs must be reduced */
#define LL_ADDMOD(r1, r0, a1, a0, b1, b0, n1, n0)                           \
    do {                                                                    \
        ulong __s2 = 0, __s1, __s0;                                         \
        add_sssaaaaaa(__s2, __s1, __s0, 0, a1, a0, 0, b1, b0);              \
        if (__s2 || __s1 > (n1) || (__s1 == (n1) && __s0 >= (n0)))          \
            sub_ddmmss(__s1, __s0, __s1, __s0, (n1), (n0));                 \
        (r1) = __s1;                                                        \
        (r0) = __s0;                                                        \
    } while (0)

/* subtract {b1,b0} from {a1,a0} modulo {n1,n0}; inputs must be reduced */
#define LL_SUBMOD(r1, r0, a1, a0, b1, b0, n1, n0)                           \
    do {                                                                    \
        ulong __d1, __d0;                                                   \
        if ((a1) > (b1) || ((a1) == (b1) && (a0) >= (b0)))                  \
            sub_ddmmss(__d1, __d0, a1, a0, b1, b0);                         \
        else                                                                \
        {                                                                   \
            sub_ddmmss(__d1, __d0, (n1), (n0), b1, b0);                     \
            add_ssaaaa(__d1, __d0, __d1, __d0, a1, a0);                     \
        }                                                                   \
        (r1) = __d1;                                                        \
        (r0) = __d0;                                                        \
    } while (0)

/*
   The Montgomery constant is -n^(-1) mod 2^FLINT_BITS; n_binvert gives
   n^(-1) mod 2^FLINT_BITS (Hurchalla's algorithm with a lookup table).
*/
#define LL_NEG_NINV(n0) (-n_binvert(n0))

#endif

#endif

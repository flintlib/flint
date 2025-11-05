/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_MACROS_DFT_H
#define N_FFT_MACROS_DFT_H

#include "longlong.h"      /* for umul_ppmm */
#include "ulong_extras.h"  /* for mulmod_shoup* functions */

/*---------*/
/* helpers */
/*---------*/

/** Shoup's modular multiplication with precomputation, lazy
 * (does not perform the excess correction step)
 *  --> computes either r or r+n and store it is res, where r = (a*b) % n
 *  --> a_pr is the precomputation for n, p_hi and p_lo are temporaries
 */
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n)             \
do {                                                          \
    ulong p_hi, p_lo;                                         \
    umul_ppmm(p_hi, p_lo, (a_pr), (b));                       \
    res = (a) * (b) - p_hi * (n);                             \
} while(0)

/*------------------*/
/* length 2, node 0 */
/*------------------*/

/** Butterfly radix 2
 * * In-place transform:                    [1  1]
 *                         [a  b] <- [a  b] [1 -1]
 * * n is the modulus, n2 is 2*n
 * * lazy_1_2:    in [0..n) / out [0..2n) / max < 2n
 * * lazy_22_24:  in [0..2n) x [0..2n) / out [0..2n) x [0..4n) / max < 4n
 * * lazy_42_44:  in [0..4n) x [0..2n) / out [0..4n) x [0..4n) / max < 4n
 */
#define DFT2_LAZY_1_2(a, b, n) \
do {                           \
    ulong tmp;                 \
    tmp = (b);                 \
    (b) = (a) + (n) - tmp;     \
    (a) = (a) + tmp;           \
} while(0)

#define DFT2_LAZY_22_24(a, b, n2) \
do {                              \
    ulong tmp;                    \
    tmp = (b);                    \
    (b) = (a) + (n2) - tmp;       \
    (a) = (a) + tmp;              \
    if ((a) >= (n2))              \
        (a) -= (n2);              \
} while(0)

#define DFT2_LAZY_42_44(a, b, n2)          \
do {                                       \
    ulong tmp;                             \
    tmp = (a);                             \
    if (tmp >= (n2))                       \
        tmp -= (n2);         /* [0..2n) */ \
    (a) = tmp + (b);         /* [0..4n) */ \
    (b) = tmp + (n2) - (b);  /* [0..4n) */ \
} while(0)

/*----------------------------------------------*/
/* length 2, general node                       */
/* (Cooley-Tukey & Gentleman-Sande butterflies) */
/*----------------------------------------------*/

/** Cooley-Tukey butterfly:
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [w -w]
 * * n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
 * * can be seen as evaluation at points w and -w of a+b*x
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define DFT2_NODE_LAZY_4_4(a, b, w, w_pr, n, n2) \
do {                                             \
    ulong u, v;                                  \
    u = (a);                                     \
    if (u >= (n2))                               \
        u -= (n2);  /* [0..2n) */                \
    v = (b);                                     \
    N_MULMOD_PRECOMP_LAZY(v, w, v, w_pr, n);     \
    (a) = u + v;                                 \
    (b) = u + (n2) - v;                          \
} while(0)

/** Gentleman-Sande butterfly:
 * * In-place transform
 *                            [1  w]
 *           [a  b] <- [a  b] [1 -w]
 * * n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
 * * can be seen as degree-1 interpolation at points iw = 1 / w and -iw, up to
 * a scaling by 1/2, since the inverse of [1  w] is  1/2 * [ 1   1]
 *                                        [1 -w]           [iw -iw]
 * * lazy_22: in [0..2n) / out [0..2n) / max < 4n
 */
#define IDFT2_NODE_LAZY_2_2(a, b, w, w_pr,        \
                            n, n2)                \
do {                                              \
    ulong tmp;                                    \
    tmp = (a) + (n2) - (b);     /* [0..4n) */     \
    (a) = (a) + (b);            /* [0..4n) */     \
    if ((a) >= (n2))                              \
        (a) -= (n2);            /* [0..2n) */     \
    N_MULMOD_PRECOMP_LAZY((b), w, tmp, w_pr, n);  \
                        /* --> (b) in [0..2n) */  \
} while(0)

/*------------------*/
/* length 4, node 0 */
/*------------------*/

/** 4-point FFT evaluation
 * * In-place transform
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
 *                              [1  0  1  0] [1  1  0  0] 
 *              == [a  b  c  d] [0  1  0  I] [1 -1  0  0] 
 *                              [1  0 -1  0] [0  0  1  1] 
 *                              [0  1  0 -I] [0  0  1 -1] 
 * * Corresponds to reducing down the tree with nodes
 *                       x^4 - 1
 *                     /         \
 *             x^2 - 1             x^2 + 1
 *             /     \             /     \
 *         x - 1     x + 1     x - I     x + I
 *  where I is typically a square root of -1
 *  (but this property is not exploited)
 * * n is the modulus, n2 is 2*n
 *   I_pr is the precomputed data for multiplication by I mod n
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 */
#define DFT4_LAZY_1_4(a, b, c, d, I, I_pr, n, n2)               \
do {                                                            \
    const ulong v0 = (a);                                       \
    const ulong v1 = (b);                                       \
    const ulong v2 = (c);                                       \
    const ulong v3 = (d);                                       \
    ulong v4 = v0 + v2;                     /* < 2*n */         \
    ulong v5 = v0 + (n) - v2;               /* < 2*n */         \
    ulong v6 = v1 + v3;                     /* < 2*n */         \
    ulong v7;                                                   \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v1 + (n) - v3, (I_pr), (n)); \
    (a) = v4 + v6;                          /* < 4*n */         \
    (b) = v4 + (n2) - v6;                   /* < 4*n */         \
    (c) = v5 + v7;                          /* < 4*n */         \
    (d) = v5 + (n2) - v7;                   /* < 4*n */         \
} while(0)

#define DFT4_LAZY_2_4(a, b, c, d, I, I_pr, n, n2)                 \
do {                                                              \
    const ulong v0 = (a);                                         \
    const ulong v1 = (b);                                         \
    const ulong v2 = (c);                                         \
    const ulong v3 = (d);                                         \
    ulong v4 = v0 + v2;                      /* < 4*n */          \
    if (v4 >= (n2))                                               \
        v4 -= (n2);                          /* < 2*n */          \
    ulong v5 = v0 + (n2) - v2;               /* < 4*n */          \
    if (v5 >= (n2))                                               \
        v5 -= (n2);                          /* < 2*n */          \
    ulong v6 = v1 + v3;                      /* < 4*n */          \
    if (v6 >= (n2))                                               \
        v6 -= (n2);                          /* < 2*n */          \
    ulong v7;                                                     \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v1 + (n2) - v3, (I_pr), (n));  \
    (a) = v4 + v6;                           /* < 4*n */          \
    (b) = v4 + (n2) - v6;                    /* < 4*n */          \
    (c) = v5 + v7;                           /* < 4*n */          \
    (d) = v5 + (n2) - v7;                    /* < 4*n */          \
} while(0)

/** 4-point FFT interpolation
 * * In-place transform
 *                              [1  1  1  1]
 *                              [1 -1  1 -1]
 * [a  b  c  d] <- [a  b  c  d] [1 -I -1  I]
 *                              [1  I -1 -I]
 *                              [1  1  0  0] [1  0  1  0]
 *              == [a  b  c  d] [1 -1  0  0] [0  1  0  1]
 *                              [0  0  1  I] [1  0 -1  0]
 *                              [0  0  1 -I] [0  1  0 -1]
 *
 * * If I**2 == -1, this matrix is the inverse of the one above; this
 * corresponds to interpolation at 1, -1, I, -I, up to scaling by 1/4; or to
 * going up the tree with nodes
 *                       x^4 - 1
 *                     /         \
 *             x^2 - 1             x^2 + 1
 *             /     \             /     \
 *         x - 1     x + 1     x - I     x + I
 * * n is the modulus, n2 is 2*n
 *   I_pr is the precomputed data for multiplication by I mod n
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_4222_4: a in [0..4n), b,c,d in [0..2n) / out [0..4n) / max < 4n
 */
#define IDFT4_LAZY_1_4(a, b, c, d, I, I_pr, n, n2)               \
do {                                                             \
    const ulong v0 = (a);                                        \
    const ulong v1 = (b);                                        \
    const ulong v2 = (c);                                        \
    const ulong v3 = (d);                                        \
    ulong v4 = v0 + v1;                         /* < 2*n */      \
    ulong v5 = v0 + (n) - v1;                   /* < 2*n */      \
    ulong v6 = v2 + v3;                         /* < 2*n */      \
    ulong v7;                                                    \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n) - v3, (I_pr), (n));  \
    (a) = v4 + v6;                              /* < 4*n */      \
    (b) = v5 + v7;                              /* < 4*n */      \
    (c) = v4 + (n2) - v6;                       /* < 4*n */      \
    (d) = v5 + (n2) - v7;                       /* < 4*n */      \
} while(0)

#define IDFT4_LAZY_4222_4(a, b, c, d, I, I_pr, n, n2)             \
do {                                                              \
    ulong v0 = (a);                                               \
    const ulong v1 = (b);                                         \
    const ulong v2 = (c);                                         \
    const ulong v3 = (d);                                         \
    if (v0 >= (n2))                                               \
        v0 -= (n2);                             /* < 2*n */       \
    ulong v4 = v0 + v1;                         /* < 4*n */       \
    if (v4 >= (n2))                                               \
        v4 -= (n2);                             /* < 2*n */       \
    ulong v5 = v0 + (n2) - v1;                  /* < 4*n */       \
    if (v5 >= (n2))                                               \
        v5 -= (n2);                             /* < 2*n */       \
    ulong v6 = v2 + v3;                         /* < 4*n */       \
    if (v6 >= (n2))                                               \
        v6 -= (n2);                             /* < 2*n */       \
    ulong v7;                                                     \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n2) - v3, (I_pr), (n));  \
    (a) = v4 + v6;                              /* < 4*n */       \
    (b) = v5 + v7;                              /* < 4*n */       \
    (c) = v4 + (n2) - v6;                       /* < 4*n */       \
    (d) = v5 + (n2) - v7;                       /* < 4*n */       \
} while(0)

/*------------------------*/
/* length 4, general node */
/*------------------------*/

/** 4-point FFT, evaluation, general node
 * * In-place transform
 *                              [ 1          1       1       1]
 *                              [w2        -w2      w3     -w3]
 * [a  b  c  d] <- [a  b  c  d] [w1         w1     -w1     -w1]
 *                              [w1*w2  -w1*w2  -w1*w3   w1*w3]
 * * Corresponds to reducing down the tree with nodes
 *                        x^4 - w1**2
 *                      /             \
 *             x^2 - w1                 x^2 + w1
 *             /      \                 /      \
 *        x - w2      x + w2       x - w3      x + w3
 * typically w2**2 == w1 and w3 == I*w2 (hence w3**2 == -w1) so that the above
 * is a Vandermonde matrix and this tree really is the subproduct tree built
 * from the four roots w2, -w2, I*w2, -I*w2 of x**4 - w1
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define DFT4_NODE_LAZY_4_4(a, b, c, d,                        \
                           w1, w1_pr, w2, w2_pr, w3, w3_pr,   \
                           n, n2)                             \
do {                                                          \
    ulong tmp;                                                \
    ulong u0 = (a);                                           \
    ulong u1 = (b);                                           \
    ulong u2 = (c);                                           \
    ulong u3 = (d);                                           \
    if (u0 >= n2)                                             \
        u0 -= n2;                                             \
    if (u1 >= n2)                                             \
        u1 -= n2;                                             \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u2, w1, u2, w1_pr, n);              \
    tmp = u0;                                                 \
    u0 = u0 + u2;                    /* [0..4n) */            \
    u2 = tmp + n2 - u2;              /* [0..4n) */            \
    if (u0 >= n2)                                             \
        u0 -= n2;                    /* [0..2n) */            \
    if (u2 >= n2)                                             \
        u2 -= n2;                    /* [0..2n) */            \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u3, w1, u3, w1_pr, n);              \
    tmp = u1;                                                 \
    u1 = u1 + u3;                    /* [0..4n) */            \
    u3 = tmp + n2 - u3;              /* [0..4n) */            \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n);              \
    (a) = u0 + u1;                   /* [0..4n) */            \
    (b) = u0 + n2 - u1;              /* [0..4n) */            \
                                                              \
    N_MULMOD_PRECOMP_LAZY(u3, w3, u3, w3_pr, n);              \
    (c) = u2 + u3;                    /* [0..4n) */           \
    (d) = u2 + n2 - u3;              /* [0..4n) */            \
} while(0)

/** 4-point FFT, interpolation, general node
 * * In-place transform
 *                              [ 1   iw2   iw1    iw1*iw2]
 *                              [ 1  -iw2   iw1   -iw1*iw2]
 * [a  b  c  d] <- [a  b  c  d] [ 1   iw3  -iw1   -iw1*iw3]
 *                              [ 1  -iw3  -iw1    iw1*iw3]
 *                              [1  iw2  0    0] [1   0   w1   0]
 *              == [a  b  c  d] [1 -iw2  0    0] [0   1    0  w1]
 *                              [0    0  1  iw3] [1   0  -w1   0]
 *                              [0    0  1 -iw3] [0   1    0 -w1]
 * * Corresponds, up to scaling by 1/4, to going up the tree with nodes
 *                        x^4 - w1**2
 *                      /             \
 *             x^2 - w1                 x^2 + w1
 *             /      \                 /      \
 *        x - w2      x + w2       x - w3      x + w3
 * typically w2**2 == w1 and w3 == I*w2 (hence w3**2 == -w1) so that the above
 * is the inverse of a Vandermonde matrix and this tree really is the
 * subproduct tree built from the four roots w2, -w2, I*w2, -I*w2 of x**4 - w1
 * * lazy_1_2: in [0..n) / out [0..2n) / max < 4n
 * * lazy_2_2: in [0..2n) / out [0..2n) / max < 4n
 */
#define IDFT4_NODE_LAZY_2_2(a, b, c, d,                              \
                            w1, w1_pr, w2, w2_pr, w3, w3_pr,         \
                            n, n2)                                   \
do {                                                                 \
    const ulong v0 = (a);                                            \
    const ulong v1 = (b);                                            \
    const ulong v2 = (c);                                            \
    const ulong v3 = (d);                                            \
    ulong v4 = v0 + v1;                       /* < 4*n */            \
    if (v4 >= (n2))                                                  \
        v4 -= (n2);                           /* < 2*n */            \
    ulong v5;                                                        \
    N_MULMOD_PRECOMP_LAZY(v5, (w2), v0 + (n2) - v1, (w2_pr), (n));   \
    ulong v6 = v2 + v3;                       /* < 4*n */            \
    if (v6 >= (n2))                                                  \
        v6 -= (n2);                           /* < 2*n */            \
    ulong v7;                                                        \
    N_MULMOD_PRECOMP_LAZY(v7, (w3), v2 + (n2) - v3, (w3_pr), (n));   \
                                                                     \
    (a) = v4 + v6;                                                   \
    if ((a) >= (n2))                                                 \
        (a) -= (n2);                           /* < 2*n */           \
    (b) = v5 + v7;                                                   \
    if ((b) >= (n2))                                                 \
        (b) -= (n2);                           /* < 2*n */           \
    N_MULMOD_PRECOMP_LAZY((c), (w1), v4 + (n2) - v6, (w1_pr), (n));  \
    N_MULMOD_PRECOMP_LAZY((d), (w1), v5 + (n2) - v7, (w1_pr), (n));  \
} while(0)

#define IDFT4_NODE_LAZY_1_2(a, b, c, d,                              \
                            w1, w1_pr, w2, w2_pr, w3, w3_pr,         \
                            n, n2)                                   \
do {                                                                 \
    const ulong v0 = (a);                                            \
    const ulong v1 = (b);                                            \
    const ulong v2 = (c);                                            \
    const ulong v3 = (d);                                            \
    ulong v4 = v0 + v1;                       /* < 2*n */            \
    ulong v5;                                                        \
    N_MULMOD_PRECOMP_LAZY(v5, (w2), v0 + (n) - v1, (w2_pr), (n));    \
    ulong v6 = v2 + v3;                       /* < 2*n */            \
    ulong v7;                                                        \
    N_MULMOD_PRECOMP_LAZY(v7, (w3), v2 + (n) - v3, (w3_pr), (n));    \
                                                                     \
    (a) = v4 + v6;                            /* < 4*n */            \
    if ((a) >= (n2))                                                 \
        (a) -= (n2);                           /* < 2*n */           \
    (b) = v5 + v7;                            /* < 4*n */            \
    if ((b) >= (n2))                                                 \
        (b) -= (n2);                           /* < 2*n */           \
    N_MULMOD_PRECOMP_LAZY((c), (w1), v4 + (n2) - v6, (w1_pr), (n));  \
    N_MULMOD_PRECOMP_LAZY((d), (w1), v5 + (n2) - v7, (w1_pr), (n));  \
} while(0)

/*------------------*/
/* length 8, node 0 */
/*------------------*/

/** 8-point FFT, evaluation
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7 into its evaluations
 *       p(1), p(-1), p(I), p(-I), p(J), p(-J), p(I*J), p(-I*J)
 * i.e. the evaluations at all 8-th roots of unity J**k for 0 <= k < 8 in
 * bit-reversed order
 * * Recall: [F->tab_w[2*k] for k in range(4)] == [1, I, J, IJ])
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 */
#define DFT8_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7, \
                      n, n2, tab_w)                   \
do {                                                  \
    DFT2_LAZY_1_2(p0, p4, n);                         \
    DFT2_LAZY_1_2(p1, p5, n);                         \
    DFT2_LAZY_1_2(p2, p6, n);                         \
    DFT2_LAZY_1_2(p3, p7, n);                         \
                                                      \
    DFT4_LAZY_2_4(p0, p1, p2, p3,                     \
                tab_w[2], tab_w[3],                   \
                n, n2);                               \
    /* could use a lazy_2_4 variant of the  */        \
    /* next one, but the gain is negligible */        \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                \
                     tab_w[2], tab_w[3],              \
                     tab_w[4], tab_w[5],              \
                     tab_w[6], tab_w[7],              \
                     n, n2);                          \
} while(0)

#define DFT8_LAZY_2_4(p0, p1, p2, p3, p4, p5, p6, p7,  \
                      n, n2, tab_w)                    \
do {                                                   \
    DFT2_LAZY_22_24(p0, p4, n2);                       \
    DFT2_LAZY_22_24(p1, p5, n2);                       \
    DFT2_LAZY_22_24(p2, p6, n2);                       \
    DFT2_LAZY_22_24(p3, p7, n2);                       \
                                                       \
    DFT4_LAZY_2_4(p0, p1, p2, p3,                      \
                tab_w[2], tab_w[3],                    \
                n, n2);                                \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                 \
                     tab_w[2], tab_w[3],               \
                     tab_w[4], tab_w[5],               \
                     tab_w[6], tab_w[7],               \
                     n, n2);                           \
} while(0)

/** 8-point FFT, interpolation
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as the evaluations
 *       [p(1), p(-1), p(I), p(-I), p(J), p(-J), p(I*J), p(-I*J)]
 * at all 8-th roots of unity J**k for 0 <= k < 8 in bit-reversed order
 * of a polynomial p(x) of degree < 8, into the coefficients of this polynomial 
 * * Recall: [F->tab_w[2*k] for k in range(4)] == [1, I, J, IJ])
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define IDFT8_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7, \
                       n, n2, tab_w)                   \
do {                                                   \
    IDFT4_LAZY_1_4(p0, p1, p2, p3,                     \
                   tab_w[2], tab_w[3],                 \
                   n, n2);                             \
    IDFT4_NODE_LAZY_1_2(p4, p5, p6, p7,                \
                tab_w[2], tab_w[3],                    \
                tab_w[4], tab_w[5],                    \
                tab_w[6], tab_w[7],                    \
                n, n2);                                \
                                                       \
    DFT2_LAZY_42_44(p0, p4, n2);                       \
    DFT2_LAZY_42_44(p1, p5, n2);                       \
    DFT2_LAZY_42_44(p2, p6, n2);                       \
    DFT2_LAZY_42_44(p3, p7, n2);                       \
} while(0)

/*------------------------*/
/* length 8, general node */
/*------------------------*/

/** 8-point FFT, evaluation, general node
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7, into its evaluations
 *       p(w0), p(-w0), p(w1), p(-w1), p(w2), p(-w2), p(w3), p(-w3)
 * where w_k = F->tab_w[8*node + 2*k] for 0 <= k < 4
 * * By construction these 8 evaluation points are the 8 roots of the
 * polynomial x**8 - F->tab_w[node]
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,    \
                           node, n, n2, tab_w)                \
do {                                                          \
    const ulong w = tab_w[2*(node)];                          \
    const ulong w_pr = tab_w[2*(node)+1];                     \
    DFT2_NODE_LAZY_4_4(p0, p4, w, w_pr, n, n2);               \
    DFT2_NODE_LAZY_4_4(p1, p5, w, w_pr, n, n2);               \
    DFT2_NODE_LAZY_4_4(p2, p6, w, w_pr, n, n2);               \
    DFT2_NODE_LAZY_4_4(p3, p7, w, w_pr, n, n2);               \
                                                              \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                        \
                       tab_w[4*(node)], tab_w[4*(node)+1],    \
                       tab_w[8*(node)], tab_w[8*(node)+1],    \
                       tab_w[8*(node)+2], tab_w[8*(node)+3],  \
                       n, n2);                                \
                                                              \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                        \
                       tab_w[4*(node)+2], tab_w[4*(node)+3],  \
                       tab_w[8*(node)+4], tab_w[8*(node)+5],  \
                       tab_w[8*(node)+6], tab_w[8*(node)+7],  \
                       n, n2);                                \
} while(0)

/** 8-point FFT, interpolation, general node
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as the evaluations
 *       [p(w0), p(-w0), p(w1), p(-w1), p(w2), p(-w2), p(w3), p(-w3)]
 * where w_k = F->tab_w[8*node + 2*k] for 0 <= k < 4 of a polynomial p(x) of
 * degree < 8, into the coefficients of this polynomial 
 * * By construction these 8 evaluation points are the 8 roots of the
 * polynomial x**8 - F->tab_w[node]
 * * lazy_1_2: in [0..n) / out [0..2n) / max < 4n
 * * lazy_2_2: in [0..2n) / out [0..2n) / max < 4n
 */
#define IDFT8_NODE_LAZY_1_2(p0, p1, p2, p3, p4, p5, p6, p7,    \
                            node, n, n2, tab_w)                \
do {                                                           \
    const ulong w = tab_w[2*(node)];                           \
    const ulong w_pr = tab_w[2*(node)+1];                      \
                                                               \
    IDFT4_NODE_LAZY_1_2(p0, p1, p2, p3,                        \
                        tab_w[4*(node)], tab_w[4*(node)+1],    \
                        tab_w[8*(node)], tab_w[8*(node)+1],    \
                        tab_w[8*(node)+2], tab_w[8*(node)+3],  \
                        n, n2);                                \
                                                               \
    IDFT4_NODE_LAZY_1_2(p4, p5, p6, p7,                        \
                        tab_w[4*(node)+2], tab_w[4*(node)+3],  \
                        tab_w[8*(node)+4], tab_w[8*(node)+5],  \
                        tab_w[8*(node)+6], tab_w[8*(node)+7],  \
                        n, n2);                                \
                                                               \
    IDFT2_NODE_LAZY_2_2(p0, p4, w, w_pr, n, n2);               \
    IDFT2_NODE_LAZY_2_2(p1, p5, w, w_pr, n, n2);               \
    IDFT2_NODE_LAZY_2_2(p2, p6, w, w_pr, n, n2);               \
    IDFT2_NODE_LAZY_2_2(p3, p7, w, w_pr, n, n2);               \
} while(0)

/*-------------------*/
/* length 16, node 0 */
/*-------------------*/

/** 16-point FFT, evaluation
 * * In-place transform p of length 16, seen as a polynomial
 * p(x) = p0 + p1*x + ... + p15*x**15, into its evaluations
 * at all 16-th roots of unity 1, -1, I, -I... (bit-reversed order)
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 */
#define DFT16_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7,           \
                       p8, p9, p10, p11, p12, p13, p14, p15,     \
                       n, n2, tab_w)                             \
do {                                                             \
    DFT4_LAZY_1_4(p0, p4, p8, p12, tab_w[2], tab_w[3], n, n2);   \
    if (p0 >= n2)                                                \
        p0 -= n2;                                                \
    DFT4_LAZY_1_4(p1, p5, p9, p13, tab_w[2], tab_w[3], n, n2);   \
    if (p1 >= n2)                                                \
        p1 -= n2;                                                \
    DFT4_LAZY_1_4(p2, p6, p10, p14, tab_w[2], tab_w[3], n, n2);  \
    if (p2 >= n2)                                                \
        p2 -= n2;                                                \
    DFT4_LAZY_1_4(p3, p7, p11, p15, tab_w[2], tab_w[3], n, n2);  \
    if (p3 >= n2)                                                \
        p3 -= n2;                                                \
                                                                 \
    /* next line requires < 2n,        */                        \
    /* hence the four reductions above */                        \
    DFT4_LAZY_2_4(p0, p1, p2, p3, tab_w[2], tab_w[3], n, n2);    \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                           \
                       tab_w[2], tab_w[3],                       \
                       tab_w[4], tab_w[5],                       \
                       tab_w[6], tab_w[7],                       \
                       n, n2);                                   \
    DFT4_NODE_LAZY_4_4(p8, p9, p10, p11,                         \
                       tab_w[4], tab_w[5],                       \
                       tab_w[8], tab_w[9],                       \
                       tab_w[10], tab_w[11],                     \
                       n, n2);                                   \
    DFT4_NODE_LAZY_4_4(p12, p13, p14, p15,                       \
                       tab_w[6], tab_w[7],                       \
                       tab_w[12], tab_w[13],                     \
                       tab_w[14], tab_w[15],                     \
                       n, n2);                                   \
} while(0)

#define DFT16_LAZY_2_4(p0, p1, p2, p3, p4, p5, p6, p7,           \
                       p8, p9, p10, p11, p12, p13, p14, p15,     \
                       n, n2, tab_w)                             \
do {                                                             \
    DFT4_LAZY_2_4(p0, p4, p8, p12, tab_w[2], tab_w[3], n, n2);   \
    if (p0 >= n2)                                                \
        p0 -= n2;                                                \
    DFT4_LAZY_2_4(p1, p5, p9, p13, tab_w[2], tab_w[3], n, n2);   \
    if (p1 >= n2)                                                \
        p1 -= n2;                                                \
    DFT4_LAZY_2_4(p2, p6, p10, p14, tab_w[2], tab_w[3], n, n2);  \
    if (p2 >= n2)                                                \
        p2 -= n2;                                                \
    DFT4_LAZY_2_4(p3, p7, p11, p15, tab_w[2], tab_w[3], n, n2);  \
    if (p3 >= n2)                                                \
        p3 -= n2;                                                \
                                                                 \
    /* next line requires < 2n,        */                        \
    /* hence the four reductions above */                        \
    DFT4_LAZY_2_4(p0, p1, p2, p3, tab_w[2], tab_w[3], n, n2);    \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                           \
                       tab_w[2], tab_w[3],                       \
                       tab_w[4], tab_w[5],                       \
                       tab_w[6], tab_w[7],                       \
                       n, n2);                                   \
    DFT4_NODE_LAZY_4_4(p8, p9, p10, p11,                         \
                       tab_w[4], tab_w[5],                       \
                       tab_w[8], tab_w[9],                       \
                       tab_w[10], tab_w[11],                     \
                       n, n2);                                   \
    DFT4_NODE_LAZY_4_4(p12, p13, p14, p15,                       \
                       tab_w[6], tab_w[7],                       \
                       tab_w[12], tab_w[13],                     \
                       tab_w[14], tab_w[15],                     \
                       n, n2);                                   \
} while(0)

/** 16-point FFT, interpolation
 * * In-place transform p of length 16, seen as the evaluations at all 16-th
 * roots of unity  1, -1, I, -I... (bit-reversed order) of a polynomial p(x) of
 * degree < 16, into the coefficients of this polynomial 
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define IDFT16_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7,              \
                        p8, p9, p10, p11, p12, p13, p14, p15,        \
                        n, n2, tab_w)                                \
do {                                                                 \
    IDFT4_LAZY_1_4(p0, p1, p2, p3, tab_w[2], tab_w[3], n, n2);       \
    IDFT4_NODE_LAZY_1_2(p4, p5, p6, p7,                              \
                        tab_w[2], tab_w[3],                          \
                        tab_w[4], tab_w[5],                          \
                        tab_w[6], tab_w[7],                          \
                        n, n2);                                      \
    IDFT4_NODE_LAZY_1_2(p8, p9, p10, p11,                            \
                        tab_w[4], tab_w[5],                          \
                        tab_w[8], tab_w[9],                          \
                        tab_w[10], tab_w[11],                        \
                        n, n2);                                      \
    IDFT4_NODE_LAZY_1_2(p12, p13, p14, p15,                          \
                        tab_w[6], tab_w[7],                          \
                        tab_w[12], tab_w[13],                        \
                        tab_w[14], tab_w[15],                        \
                        n, n2);                                      \
                                                                     \
    IDFT4_LAZY_4222_4(p0, p4, p8, p12, tab_w[2], tab_w[3], n, n2);   \
    IDFT4_LAZY_4222_4(p1, p5, p9, p13, tab_w[2], tab_w[3], n, n2);   \
    IDFT4_LAZY_4222_4(p2, p6, p10, p14, tab_w[2], tab_w[3], n, n2);  \
    IDFT4_LAZY_4222_4(p3, p7, p11, p15, tab_w[2], tab_w[3], n, n2);  \
} while(0)

/*-------------------------*/
/* length 16, general node */
/*-------------------------*/

/** 16-point FFT, evaluation, general node
 * * In-place transform p of length 16, seen as a polynomial
 * p(x) = p0 + p1*x + ... + p15*x**15, into its evaluations at 
 *       p(w0), p(-w0), p(w1), p(-w1), ..., p(w7), p(-w7)
 * where w_k = F->tab_w[16*node + 2*k] for 0 <= k < 8
 * * By construction these 16 evaluation points are the 16 roots of the
 * polynomial x**16 - F->tab_w[node]
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define DFT16_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,       \
                            p8, p9, p10, p11, p12, p13, p14, p15, \
                            node, n, n2, tab_w)                   \
do {                                                              \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                          \
                                                                  \
    w2 = tab_w[2*node];                                           \
    w2pre = tab_w[2*node+1];                                      \
    w = tab_w[4*node];                                            \
    wpre = tab_w[4*node+1];                                       \
    Iw = tab_w[4*node+2];                                         \
    Iwpre = tab_w[4*node+3];                                      \
                                                                  \
    DFT4_NODE_LAZY_4_4(p0, p4, p8, p12,                           \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
    DFT4_NODE_LAZY_4_4(p1, p5, p9, p13,                           \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
    DFT4_NODE_LAZY_4_4(p2, p6, p10, p14,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
    DFT4_NODE_LAZY_4_4(p3, p7, p11, p15,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
                                                                  \
    w2 = tab_w[8*node];                                           \
    w2pre = tab_w[8*node+1];                                      \
    w = tab_w[16*node];                                           \
    wpre = tab_w[16*node+1];                                      \
    Iw = tab_w[16*node+2];                                        \
    Iwpre = tab_w[16*node+3];                                     \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                            \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
                                                                  \
    w2 = tab_w[8*node+2];                                         \
    w2pre = tab_w[8*node+3];                                      \
    w = tab_w[16*node+4];                                         \
    wpre = tab_w[16*node+5];                                      \
    Iw = tab_w[16*node+6];                                        \
    Iwpre = tab_w[16*node+7];                                     \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                            \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
                                                                  \
    w2 = tab_w[8*node+4];                                         \
    w2pre = tab_w[8*node+5];                                      \
    w = tab_w[16*node+8];                                         \
    wpre = tab_w[16*node+9];                                      \
    Iw = tab_w[16*node+10];                                       \
    Iwpre = tab_w[16*node+11];                                    \
    DFT4_NODE_LAZY_4_4(p8, p9, p10, p11,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
                                                                  \
    w2 = tab_w[8*node+6];                                         \
    w2pre = tab_w[8*node+7];                                      \
    w = tab_w[16*node+12];                                        \
    wpre = tab_w[16*node+13];                                     \
    Iw = tab_w[16*node+14];                                       \
    Iwpre = tab_w[16*node+15];                                    \
    DFT4_NODE_LAZY_4_4(p12, p13, p14, p15,                        \
                       w2, w2pre, w, wpre, Iw, Iwpre,             \
                       n, n2);                                    \
} while(0)

/** 16-point FFT, interpolation, general node
 * * In-place transform p of length 16, seen as the evaluations at 
 *       w0, -w0, w1, -w1, ..., w7, -w7
 * where w_k = F->tab_w[16*node + 2*k] for 0 <= k < 8
 * of a polynomial of degree < 16, into the coefficients of this polynomial 
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define IDFT16_NODE_LAZY_1_2(p0, p1, p2, p3, p4, p5, p6, p7,       \
                             p8, p9, p10, p11, p12, p13, p14, p15, \
                             node, n, n2, tab_w)                   \
do {                                                               \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                           \
                                                                   \
    w2 = tab_w[8*node];                                            \
    w2pre = tab_w[8*node+1];                                       \
    w = tab_w[16*node];                                            \
    wpre = tab_w[16*node+1];                                       \
    Iw = tab_w[16*node+2];                                         \
    Iwpre = tab_w[16*node+3];                                      \
    IDFT4_NODE_LAZY_1_2(p0, p1, p2, p3,                            \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
                                                                   \
    w2 = tab_w[8*node+2];                                          \
    w2pre = tab_w[8*node+3];                                       \
    w = tab_w[16*node+4];                                          \
    wpre = tab_w[16*node+5];                                       \
    Iw = tab_w[16*node+6];                                         \
    Iwpre = tab_w[16*node+7];                                      \
    IDFT4_NODE_LAZY_1_2(p4, p5, p6, p7,                            \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
                                                                   \
    w2 = tab_w[8*node+4];                                          \
    w2pre = tab_w[8*node+5];                                       \
    w = tab_w[16*node+8];                                          \
    wpre = tab_w[16*node+9];                                       \
    Iw = tab_w[16*node+10];                                        \
    Iwpre = tab_w[16*node+11];                                     \
    IDFT4_NODE_LAZY_1_2(p8, p9, p10, p11,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
                                                                   \
    w2 = tab_w[8*node+6];                                          \
    w2pre = tab_w[8*node+7];                                       \
    w = tab_w[16*node+12];                                         \
    wpre = tab_w[16*node+13];                                      \
    Iw = tab_w[16*node+14];                                        \
    Iwpre = tab_w[16*node+15];                                     \
    IDFT4_NODE_LAZY_1_2(p12, p13, p14, p15,                        \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
                                                                   \
    w2 = tab_w[2*node];                                            \
    w2pre = tab_w[2*node+1];                                       \
    w = tab_w[4*node];                                             \
    wpre = tab_w[4*node+1];                                        \
    Iw = tab_w[4*node+2];                                          \
    Iwpre = tab_w[4*node+3];                                       \
                                                                   \
    IDFT4_NODE_LAZY_2_2(p0, p4, p8, p12,                           \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
    IDFT4_NODE_LAZY_2_2(p1, p5, p9, p13,                           \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
    IDFT4_NODE_LAZY_2_2(p2, p6, p10, p14,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
    IDFT4_NODE_LAZY_2_2(p3, p7, p11, p15,                          \
                       w2, w2pre, w, wpre, Iw, Iwpre,              \
                       n, n2);                                     \
} while(0)


/*-------------------*/
/* length 32, node 0 */
/*-------------------*/

/** 32-point FFT, evaluation
 * * In-place transform p of length 32, seen as a polynomial
 * p(x) = p0 + p1*x + ... + p31*x**31, into its evaluations
 * at all 32-th roots of unity 1, -1, I, -I... (bit-reversed order)
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 */
#define DFT32_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7,                           \
                       p8, p9, p10, p11, p12, p13, p14, p15,                     \
                       p16, p17, p18, p19, p20, p21, p22, p23,                   \
                       p24, p25, p26, p27, p28, p29, p30, p31,                   \
                       n, n2, tab_w)                                             \
do {                                                                             \
    DFT4_LAZY_1_4(p0, p8, p16, p24, tab_w[2], tab_w[3], n, n2);                  \
    if (p0 >= n2)                                                                \
        p0 -= n2;                                                                \
    DFT4_LAZY_1_4(p1, p9, p17, p25, tab_w[2], tab_w[3], n, n2);                  \
    if (p1 >= n2)                                                                \
        p1 -= n2;                                                                \
    DFT4_LAZY_1_4(p2, p10, p18, p26, tab_w[2], tab_w[3], n, n2);                 \
    if (p2 >= n2)                                                                \
        p2 -= n2;                                                                \
    DFT4_LAZY_1_4(p3, p11, p19, p27, tab_w[2], tab_w[3], n, n2);                 \
    if (p3 >= n2)                                                                \
        p3 -= n2;                                                                \
    DFT4_LAZY_1_4(p4, p12, p20, p28, tab_w[2], tab_w[3], n, n2);                 \
    if (p4 >= n2)                                                                \
        p4 -= n2;                                                                \
    DFT4_LAZY_1_4(p5, p13, p21, p29, tab_w[2], tab_w[3], n, n2);                 \
    if (p5 >= n2)                                                                \
        p5 -= n2;                                                                \
    DFT4_LAZY_1_4(p6, p14, p22, p30, tab_w[2], tab_w[3], n, n2);                 \
    if (p6 >= n2)                                                                \
        p6 -= n2;                                                                \
    DFT4_LAZY_1_4(p7, p15, p23, p31, tab_w[2], tab_w[3], n, n2);                 \
    if (p7 >= n2)                                                                \
        p7 -= n2;                                                                \
                                                                                 \
    /* next line requires < 2n, hence the 8 reductions above */                  \
    DFT8_LAZY_2_4(p0, p1, p2, p3, p4, p5, p6, p7, n, n2, tab_w);                 \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 1, n, n2, tab_w);   \
    DFT8_NODE_LAZY_4_4(p16, p17, p18, p19, p20, p21, p22, p23, 2, n, n2, tab_w); \
    DFT8_NODE_LAZY_4_4(p24, p25, p26, p27, p28, p29, p30, p31, 3, n, n2, tab_w); \
} while(0)

#define DFT32_LAZY_2_4(p0, p1, p2, p3, p4, p5, p6, p7,                           \
                       p8, p9, p10, p11, p12, p13, p14, p15,                     \
                       p16, p17, p18, p19, p20, p21, p22, p23,                   \
                       p24, p25, p26, p27, p28, p29, p30, p31,                   \
                       n, n2, tab_w)                                             \
do {                                                                             \
    DFT4_LAZY_2_4(p0, p8, p16, p24, tab_w[2], tab_w[3], n, n2);                  \
    if (p0 >= n2)                                                                \
        p0 -= n2;                                                                \
    DFT4_LAZY_2_4(p1, p9, p17, p25, tab_w[2], tab_w[3], n, n2);                  \
    if (p1 >= n2)                                                                \
        p1 -= n2;                                                                \
    DFT4_LAZY_2_4(p2, p10, p18, p26, tab_w[2], tab_w[3], n, n2);                 \
    if (p2 >= n2)                                                                \
        p2 -= n2;                                                                \
    DFT4_LAZY_2_4(p3, p11, p19, p27, tab_w[2], tab_w[3], n, n2);                 \
    if (p3 >= n2)                                                                \
        p3 -= n2;                                                                \
    DFT4_LAZY_2_4(p4, p12, p20, p28, tab_w[2], tab_w[3], n, n2);                 \
    if (p4 >= n2)                                                                \
        p4 -= n2;                                                                \
    DFT4_LAZY_2_4(p5, p13, p21, p29, tab_w[2], tab_w[3], n, n2);                 \
    if (p5 >= n2)                                                                \
        p5 -= n2;                                                                \
    DFT4_LAZY_2_4(p6, p14, p22, p30, tab_w[2], tab_w[3], n, n2);                 \
    if (p6 >= n2)                                                                \
        p6 -= n2;                                                                \
    DFT4_LAZY_2_4(p7, p15, p23, p31, tab_w[2], tab_w[3], n, n2);                 \
    if (p7 >= n2)                                                                \
        p7 -= n2;                                                                \
                                                                                 \
    /* next line requires < 2n, hence the 8 reductions above */                  \
    DFT8_LAZY_2_4(p0, p1, p2, p3, p4, p5, p6, p7, n, n2, tab_w);                 \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 1, n, n2, tab_w);   \
    DFT8_NODE_LAZY_4_4(p16, p17, p18, p19, p20, p21, p22, p23, 2, n, n2, tab_w); \
    DFT8_NODE_LAZY_4_4(p24, p25, p26, p27, p28, p29, p30, p31, 3, n, n2, tab_w); \
} while(0)

/** 32-point FFT, interpolation
 * * In-place transform p of length 32, seen as the evaluations at all 32-th
 * roots of unity  1, -1, I, -I... (bit-reversed order) of a polynomial p(x) of
 * degree < 32, into the coefficients of this polynomial 
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define IDFT32_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7,                           \
                        p8, p9, p10, p11, p12, p13, p14, p15,                     \
                        p16, p17, p18, p19, p20, p21, p22, p23,                   \
                        p24, p25, p26, p27, p28, p29, p30, p31,                   \
                        n, n2, tab_w)                                             \
do {                                                                              \
    IDFT8_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7, n, n2, tab_w);                 \
    IDFT8_NODE_LAZY_1_2(p8, p9, p10, p11, p12, p13, p14, p15, 1, n, n2, tab_w);   \
    IDFT8_NODE_LAZY_1_2(p16, p17, p18, p19, p20, p21, p22, p23, 2, n, n2, tab_w); \
    IDFT8_NODE_LAZY_1_2(p24, p25, p26, p27, p28, p29, p30, p31, 3, n, n2, tab_w); \
                                                                                  \
    IDFT4_LAZY_4222_4(p0, p8, p16, p24, tab_w[2], tab_w[3], n, n2);               \
    IDFT4_LAZY_4222_4(p1, p9, p17, p25, tab_w[2], tab_w[3], n, n2);               \
    IDFT4_LAZY_4222_4(p2, p10, p18, p26, tab_w[2], tab_w[3], n, n2);              \
    IDFT4_LAZY_4222_4(p3, p11, p19, p27, tab_w[2], tab_w[3], n, n2);              \
    IDFT4_LAZY_4222_4(p4, p12, p20, p28, tab_w[2], tab_w[3], n, n2);              \
    IDFT4_LAZY_4222_4(p5, p13, p21, p29, tab_w[2], tab_w[3], n, n2);              \
    IDFT4_LAZY_4222_4(p6, p14, p22, p30, tab_w[2], tab_w[3], n, n2);              \
    IDFT4_LAZY_4222_4(p7, p15, p23, p31, tab_w[2], tab_w[3], n, n2);              \
} while(0)

/*-------------------------*/
/* length 32, general node */
/*-------------------------*/

/** 32-point FFT, evaluation, general node
 * * In-place transform p of length 32, seen as a polynomial
 * p(x) = p0 + p1*x + ... + p31*x**31, into its evaluations at 
 *       p(w0), p(-w0), p(w1), p(-w1), ..., p(w15), p(-w15)
 * where w_k = F->tab_w[32*node + 2*k] for 0 <= k < 16
 * * By construction these 32 evaluation points are the 32 roots of the
 * polynomial x**32 - F->tab_w[node]
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define DFT32_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                              \
                            p8, p9, p10, p11, p12, p13, p14, p15,                        \
                            p16, p17, p18, p19, p20, p21, p22, p23,                      \
                            p24, p25, p26, p27, p28, p29, p30, p31,                      \
                            node, n, n2, tab_w)                                          \
do {                                                                                     \
    ulong w2 = tab_w[2*node];                                                            \
    ulong w2pre = tab_w[2*node+1];                                                       \
    ulong w = tab_w[4*node];                                                             \
    ulong wpre = tab_w[4*node+1];                                                        \
    ulong Iw = tab_w[4*node+2];                                                          \
    ulong Iwpre = tab_w[4*node+3];                                                       \
    DFT4_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);          \
    DFT4_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);          \
    DFT4_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    DFT4_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    DFT4_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    DFT4_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    DFT4_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    DFT4_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
                                                                                         \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);            \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);    \
    DFT8_NODE_LAZY_4_4(p16, p17, p18, p19, p20, p21, p22, p23, 4*node+2, n, n2, tab_w);  \
    DFT8_NODE_LAZY_4_4(p24, p25, p26, p27, p28, p29, p30, p31, 4*node+3, n, n2, tab_w);  \
} while(0)

/** 32-point FFT, interpolation, general node
 * * In-place transform p of length 32, seen as the evaluations at 
 *       w0, -w0, w1, -w1, ..., w15, -w15
 * where w_k = F->tab_w[32*node + 2*k] for 0 <= k < 16 of a polynomial of
 * degree < 32, into the coefficients of this polynomial 
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define IDFT32_NODE_LAZY_1_2(p0, p1, p2, p3, p4, p5, p6, p7,                              \
                             p8, p9, p10, p11, p12, p13, p14, p15,                        \
                             p16, p17, p18, p19, p20, p21, p22, p23,                      \
                             p24, p25, p26, p27, p28, p29, p30, p31,                      \
                             node, n, n2, tab_w)                                          \
do {                                                                                      \
    IDFT8_NODE_LAZY_1_2(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);            \
    IDFT8_NODE_LAZY_1_2(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);    \
    IDFT8_NODE_LAZY_1_2(p16, p17, p18, p19, p20, p21, p22, p23, 4*node+2, n, n2, tab_w);  \
    IDFT8_NODE_LAZY_1_2(p24, p25, p26, p27, p28, p29, p30, p31, 4*node+3, n, n2, tab_w);  \
                                                                                          \
    ulong w2 = tab_w[2*node];                                                             \
    ulong w2pre = tab_w[2*node+1];                                                        \
    ulong w = tab_w[4*node];                                                              \
    ulong wpre = tab_w[4*node+1];                                                         \
    ulong Iw = tab_w[4*node+2];                                                           \
    ulong Iwpre = tab_w[4*node+3];                                                        \
    IDFT4_NODE_LAZY_2_2(p0, p8, p16, p24, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);          \
    IDFT4_NODE_LAZY_2_2(p1, p9, p17, p25, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);          \
    IDFT4_NODE_LAZY_2_2(p2, p10, p18, p26, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    IDFT4_NODE_LAZY_2_2(p3, p11, p19, p27, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    IDFT4_NODE_LAZY_2_2(p4, p12, p20, p28, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    IDFT4_NODE_LAZY_2_2(p5, p13, p21, p29, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    IDFT4_NODE_LAZY_2_2(p6, p14, p22, p30, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
    IDFT4_NODE_LAZY_2_2(p7, p15, p23, p31, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);         \
} while(0)

#endif  /* N_FFT_MACROS_DFT_H */

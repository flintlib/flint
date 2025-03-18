/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_MACROS_H
#define N_FFT_MACROS_H

/*---------*/
/* helpers */
/*---------*/

/** Shoup's modular multiplication with precomputation, lazy
 * (does not perform the excess correction step)
 *  --> computes either r or r+n and store it is res, where r = (a*b) % n
 *  --> a_pr is the precomputation for n, p_hi and p_lo are temporaries
 */
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n, p_hi, p_lo) \
    do {                                                      \
        umul_ppmm(p_hi, p_lo, (a_pr), (b));                   \
        res = (a) * (b) - p_hi * (n);                         \
    } while(0)

/*---------------------*/
/* radix-2 butterflies */
/*---------------------*/

/** Butterfly radix 2
 * * in [0..n) x [0..n) / out [0..2n) x [0..2n) / max < 2n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n is the modulus, tmp is a temporary
 */
#define BUTTERFLY_LAZY12(a, b, n, tmp) \
    do {                                \
        tmp = (b);                      \
        (b) = (a) + (n) - tmp;          \
        (a) = (a) + tmp;                \
    } while(0)

/** Butterfly radix 2
 * * in [0..2n) x [0..2n) / out [0..2n) x [0..4n) / max < 4n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n2 is 2*n, tmp is a temporary
 */
#define BUTTERFLY_LAZY24(a, b, n2, tmp) \
    do {                               \
        tmp = (b);                     \
        (b) = (a) + (n2) - tmp;        \
        (a) = (a) + tmp;               \
        if ((a) >= (n2))               \
            (a) -= (n2);               \
    } while(0)

/*---------------------*/
/* radix-4 butterflies */
/*---------------------*/

/** 4-point butterfly, evaluation
 * * in [0..n) / out [0..4n) / max < 4n
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
 * * n is the modulus and n2 == 2*n, p_hi, p_lo are temporaries
 */
#define DFT4_NODE0_LAZY14(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)    \
    do {                                                             \
        const ulong v0 = (a);                                        \
        const ulong v1 = (b);                                        \
        const ulong v2 = (c);                                        \
        const ulong v3 = (d);                                        \
        ulong v4 = v0 + v2;                         /* < 2*n */      \
        ulong v5 = v0 + (n) - v2;                   /* < 2*n */      \
        ulong v6 = v1 + v3;                         /* < 2*n */      \
        ulong v7;                                                    \
        N_MULMOD_PRECOMP_LAZY(v7, (I), v1 + (n) - v3, (I_pr), (n),   \
                              p_hi, p_lo);                           \
        (a) = v4 + v6;                              /* < 4*n */      \
        (b) = v4 + (n2) - v6;                       /* < 4*n */      \
        (c) = v5 + v7;                              /* < 3*n */      \
        (d) = v5 + (n2) - v7;                       /* < 4*n */      \
    } while(0)

/** 4-point butterfly, evaluation
 * * in [0..2n) / out [0..4n) / max < 4n
 * * other than this, same specification as DFT4_NODE0_LAZY14
 */
#define DFT4_NODE0_LAZY24(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)    \
    do {                                                             \
        const ulong v0 = (a);                                        \
        const ulong v1 = (b);                                        \
        const ulong v2 = (c);                                        \
        const ulong v3 = (d);                                        \
        ulong v4 = v0 + v2;                         /* < 4*n */      \
        if (v4 >= (n2))                                              \
            v4 -= (n2);                             /* < 2*n */      \
        ulong v5 = v0 + (n2) - v2;                  /* < 4*n */      \
        if (v5 >= (n2))                                              \
            v5 -= (n2);                             /* < 2*n */      \
        ulong v6 = v1 + v3;                         /* < 4*n */      \
        if (v6 >= (n2))                                              \
            v6 -= (n2);                             /* < 2*n */      \
        ulong v7;                                                    \
        N_MULMOD_PRECOMP_LAZY(v7, (I), v1 + (n2) - v3, (I_pr), (n),  \
                              p_hi, p_lo);                           \
        (a) = v4 + v6;                              /* < 4*n */      \
        (b) = v4 + (n2) - v6;                       /* < 4*n */      \
        (c) = v5 + v7;                              /* < 4*n */      \
        (d) = v5 + (n2) - v7;                       /* < 4*n */      \
    } while(0)


/** 4-point butterfly, interpolation
 * * in [0..n) / out [0..4n) / max < 4n
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
 */
#define IDFT4_NODE0_LAZY12(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo) \
do {                                                               \
    const ulong v0 = (a);                                          \
    const ulong v1 = (b);                                          \
    const ulong v2 = (c);                                          \
    const ulong v3 = (d);                                          \
    ulong v4 = v0 + v1;                         /* < 2*n */        \
    ulong v5 = v0 + (n) - v1;                   /* < 2*n */        \
    ulong v6 = v2 + v3;                         /* < 2*n */        \
    ulong v7;                                                      \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n) - v3, (I_pr), (n),     \
                            p_hi, p_lo);        /* < 2*n */        \
    (a) = v4 + v6;                              /* < 4*n */        \
    if ((a) >= (n2))                                               \
        (a) -= (n2);                            /* < 2*n */        \
    (b) = v5 + v7;                              /* < 4*n */        \
    if ((b) >= (n2))                                               \
        (b) -= (n2);                            /* < 2*n */        \
    (c) = v4 + (n2) - v6;                       /* < 4*n */        \
    if ((c) >= (n2))                                               \
        (c) -= (n2);                            /* < 2*n */        \
    (d) = v5 + (n2) - v7;                       /* < 4*n */        \
    if ((d) >= (n2))                                               \
        (d) -= (n2);                            /* < 2*n */        \
} while(0)

#define IDFT4_NODE0_LAZY14(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo) \
do {                                                               \
    const ulong v0 = (a);                                          \
    const ulong v1 = (b);                                          \
    const ulong v2 = (c);                                          \
    const ulong v3 = (d);                                          \
    ulong v4 = v0 + v1;                         /* < 2*n */        \
    ulong v5 = v0 + (n) - v1;                   /* < 2*n */        \
    ulong v6 = v2 + v3;                         /* < 2*n */        \
    ulong v7;                                                      \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n) - v3, (I_pr), (n),     \
                            p_hi, p_lo);        /* < 2*n */        \
    (a) = v4 + v6;                              /* < 4*n */        \
    (b) = v5 + v7;                              /* < 4*n */        \
    (c) = v4 + (n2) - v6;                       /* < 4*n */        \
    (d) = v5 + (n2) - v7;                       /* < 4*n */        \
} while(0)

/** 4-point butterfly, interpolation
 * * in [0..2n) / out [0..4n) / max < 4n
 * * other than this, same specification as IDFT4_NODE0_LAZY14
 */
#define IDFT4_NODE0_LAZY24(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)  \
do {                                                                \
    const ulong v0 = (a);                                           \
    const ulong v1 = (b);                                           \
    const ulong v2 = (c);                                           \
    const ulong v3 = (d);                                           \
    ulong v4 = v0 + v1;                         /* < 4*n */         \
    if (v4 >= (n2))                                                 \
        v4 -= (n2);                             /* < 2*n */         \
    ulong v5 = v0 + (n2) - v1;                  /* < 4*n */         \
    if (v5 >= (n2))                                                 \
        v5 -= (n2);                             /* < 2*n */         \
    ulong v6 = v2 + v3;                         /* < 4*n */         \
    if (v6 >= (n2))                                                 \
        v6 -= (n2);                             /* < 2*n */         \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n2) - v3, (I_pr), (n),     \
                            p_hi, p_lo);        /* < 2*n */         \
    (a) = v4 + v6;                              /* < 4*n */         \
    (b) = v5 + v7;                              /* < 4*n */         \
    (c) = v4 + (n2) - v6;                       /* < 4*n */         \
    (d) = v5 + (n2) - v7;                       /* < 4*n */         \
} while(0)

#define IDFT4_NODE0_LAZY4222(a, b, c, d, I, I_pr, n, n2, p_hi, p_lo)  \
do {                                                                \
    ulong v0 = (a);                                           \
    const ulong v1 = (b);                                           \
    const ulong v2 = (c);                                           \
    const ulong v3 = (d);                                           \
    if (v0 >= (n2))                                                 \
        v0 -= (n2);                             /* < 2*n */         \
    ulong v4 = v0 + v1;                         /* < 4*n */         \
    if (v4 >= (n2))                                                 \
        v4 -= (n2);                             /* < 2*n */         \
    ulong v5 = v0 + (n2) - v1;                  /* < 4*n */         \
    if (v5 >= (n2))                                                 \
        v5 -= (n2);                             /* < 2*n */         \
    ulong v6 = v2 + v3;                         /* < 4*n */         \
    if (v6 >= (n2))                                                 \
        v6 -= (n2);                             /* < 2*n */         \
    ulong v7;                                                       \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2 + (n2) - v3, (I_pr), (n),     \
                            p_hi, p_lo);        /* < 2*n */         \
    (a) = v4 + v6;                              /* < 4*n */         \
    (b) = v5 + v7;                              /* < 4*n */         \
    (c) = v4 + (n2) - v6;                       /* < 4*n */         \
    (d) = v5 + (n2) - v7;                       /* < 4*n */         \
} while(0)



#endif  /* N_FFT_MACROS_H */

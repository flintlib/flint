/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_MACROS_TFT_H
#define N_FFT_MACROS_TFT_H

#include "n_fft/impl_macros_dft.h"

/*----------------------------------------------*/
/* length 2, general node                       */
/* (Cooley-Tukey & Gentleman-Sande butterflies) */
/*----------------------------------------------*/

/** Cooley-Tukey butterfly, truncated:
 * * In-place transform
 *                         [1]
 *           [a] <- [a  b] [w]
 * * n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
 * * can be seen as evaluation at points w of a+b*x
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define TFT2_1_NODE_LAZY_4_4(a, b, w, w_pr, n, n2) \
do {                                               \
    ulong u, v;                                    \
    u = (a);                                       \
    if (u >= (n2))                                 \
        u -= (n2);  /* [0..2n) */                  \
    v = (b);                                       \
    N_MULMOD_PRECOMP_LAZY(v, w, v, w_pr, n);       \
    (a) = u + v;                                   \
} while(0)

/*------------------------*/
/* length 4, general node */
/*------------------------*/

/** 4-point TFT, evaluation, general node, truncated at 2
 * * same as DFT4_NODE_LAZY_4_4 but only computes first 2 evaluations in a,b;
 * the elements b,c may be modified
 * * note: TFT4_1 would be almost the same as TFT4_2,
 *     and TFT4_3 would be almost the same as DFT4
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define TFT4_2_NODE_LAZY_4_4(a, b, c, d,             \
                             w1, w1_pr, w2, w2_pr,   \
                             n, n2)                  \
do {                                                 \
    ulong u0 = (a);                                  \
    ulong u1 = (b);                                  \
    ulong u2 = (c);                                  \
    ulong u3 = (d);                                  \
    if (u0 >= n2)                                    \
        u0 -= n2;                                    \
    if (u1 >= n2)                                    \
        u1 -= n2;                                    \
                                                     \
    N_MULMOD_PRECOMP_LAZY(u2, w1, u2, w1_pr, n);     \
    u0 = u0 + u2;                    /* [0..4n) */   \
    if (u0 >= n2)                                    \
        u0 -= n2;                    /* [0..2n) */   \
                                                     \
    N_MULMOD_PRECOMP_LAZY(u3, w1, u3, w1_pr, n);     \
    u1 = u1 + u3;                    /* [0..4n) */   \
                                                     \
    N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n);     \
    (a) = u0 + u1;                   /* [0..4n) */   \
    (b) = u0 + n2 - u1;              /* [0..4n) */   \
} while(0)

/*------------------------*/
/* length 8, general node */
/*------------------------*/

/** 8-point FFT, evaluation, general node, truncated at 4
 * * same as DFT8_NODE_LAZY_4_4 but only computes first 4 evaluations in
 * p0,p1,p2,p3; other 4 elements may be modified
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define TFT8_4_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,    \
                             node, n, n2, tab_w)                \
do {                                                            \
    const ulong w = tab_w[2*(node)];                            \
    const ulong w_pr = tab_w[2*(node)+1];                       \
    TFT2_1_NODE_LAZY_4_4(p0, p4, w, w_pr, n, n2);               \
    TFT2_1_NODE_LAZY_4_4(p1, p5, w, w_pr, n, n2);               \
    TFT2_1_NODE_LAZY_4_4(p2, p6, w, w_pr, n, n2);               \
    TFT2_1_NODE_LAZY_4_4(p3, p7, w, w_pr, n, n2);               \
                                                                \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                          \
                       tab_w[4*(node)], tab_w[4*(node)+1],      \
                       tab_w[8*(node)], tab_w[8*(node)+1],      \
                       tab_w[8*(node)+2], tab_w[8*(node)+3],    \
                       n, n2);                                  \
} while(0)

/*-------------------*/
/* length 16, node 0 */
/*-------------------*/


/* TODO comment, and place where appropriate */
#define ITFT4_3_LAZY_422_2(a, b, c, d, I, I_pr, n, n2)            \
do {                                                              \
    ulong v0 = (a);                                               \
    ulong v1 = (b);                                               \
    ulong v2 = (c);                                               \
    if (v0 >= (n2))                                               \
        v0 -= (n2);                             /* < 2*n */       \
    ulong v4 = v0 + v1;                         /* < 4*n */       \
    if (v4 >= (n2))                                               \
        v4 -= (n2);                             /* < 2*n */       \
    ulong v5 = v0 + (n2) - v1;                  /* < 4*n */       \
    if (v5 >= (n2))                                               \
        v5 -= (n2);                             /* < 2*n */       \
    ulong v7;                                                     \
    N_MULMOD_PRECOMP_LAZY(v7, (I), v2, (I_pr), (n));              \
    v0 = v4 + v2;                                                 \
    v1 = v5 + v7;                                                 \
    v4 = v4 + (n2) - v2;                                          \
    v5 = v5 + (n2) - v7;                                          \
    if (v0 >= (n2))                                               \
        v0 -= (n2);                                               \
    if (v1 >= (n2))                                               \
        v1 -= (n2);                                               \
    if (v4 >= (n2))                                               \
        v4 -= (n2);                                               \
    if (v5 >= (n2))                                               \
        v5 -= (n2);                                               \
    (a) = v0;                              /* < 2*n */            \
    (b) = v1;                              /* < 2*n */            \
    (c) = v4;                              /* < 2*n */            \
    (d) = v5;                              /* < 2*n */            \
} while(0)

/** 16-point FFT, interpolation, truncated at 12
 * * same as IDFT16_LAZY_1_4 but only considers first 12 evaluations,
 * other 4 elements p12,p13,p14,p15 may be modified
 * * typically I and I_pr are the forward roots; tab_w the inverse roots
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
#define ITFT16_12_LAZY_1_4(p0, p1, p2, p3, p4, p5, p6, p7,              \
                           p8, p9, p10, p11, p12, p13, p14, p15,        \
                           I, I_pr, n, n2, tab_w)                       \
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
                                                                     \
    ITFT4_3_LAZY_422_2(p0, p4, p8, p12, tab_w[2], tab_w[3], n, n2);   \
    ITFT4_3_LAZY_422_2(p1, p5, p9, p13, tab_w[2], tab_w[3], n, n2);   \
    ITFT4_3_LAZY_422_2(p2, p6, p10, p14, tab_w[2], tab_w[3], n, n2);  \
    ITFT4_3_LAZY_422_2(p3, p7, p11, p15, tab_w[2], tab_w[3], n, n2);  \
                                                                      \
    /* reduce by (x**8 - 1) * (x**4 - I): */ \
    /* (p15*I + p11)*x^11 + (p14*I + p10)*x^10 + (p13*I + p9)*x^9 + (p12*I + p8)*x^8 */ \
    /* + (p7 + p15)*x^7 + (p6 + p14)*x^6 + (p5 + p13)*x^5 + (p4 + p12)*x^4           */ \
    /* + (-p15*I + p3)*x^3 + (-p14*I + p2)*x^2 + (-p13*I + p1)*x - p12*I + p0        */ \
    p4 += p12;                                                      \
    p5 += p13;                                                      \
    p6 += p14;                                                      \
    p7 += p15;                                                      \
    N_MULMOD_PRECOMP_LAZY(p12, (I), p12, (I_pr), (n));              \
    N_MULMOD_PRECOMP_LAZY(p13, (I), p13, (I_pr), (n));              \
    N_MULMOD_PRECOMP_LAZY(p14, (I), p14, (I_pr), (n));              \
    N_MULMOD_PRECOMP_LAZY(p15, (I), p15, (I_pr), (n));              \
    p0 = p0 + (n2) - p12;                                           \
    p1 = p1 + (n2) - p13;                                           \
    p2 = p2 + (n2) - p14;                                           \
    p3 = p3 + (n2) - p15;                                           \
    p8 += p12;                                                      \
    p9 += p13;                                                      \
    p10 += p14;                                                     \
    p11 += p15;                                                     \
} while(0)

/*-------------------------*/
/* length 16, general node */
/*-------------------------*/

/** 16-point FFT, evaluation, general node, truncated at {4,8,12}
 * * same as DFT16_NODE_LAZY_4_4 but only computes first k evaluations (for k
 * in {4, 8, 12}, 3 variants available) in the first k elements of p;
 * other 16 - k elements may be modified
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define TFT16_4_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,       \
                              p8, p9, p10, p11, p12, p13, p14, p15, \
                              node, n, n2, tab_w)                   \
do {                                                                \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                            \
                                                                    \
    w2 = tab_w[2*node];                                             \
    w2pre = tab_w[2*node+1];                                        \
    w = tab_w[4*node];                                              \
    wpre = tab_w[4*node+1];                                         \
                                                                    \
    TFT4_2_NODE_LAZY_4_4(p0, p4, p8, p12,                           \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p1, p5, p9, p13,                           \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p2, p6, p10, p14,                          \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p3, p7, p11, p15,                          \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
                                                                    \
    w2 = tab_w[8*node];                                             \
    w2pre = tab_w[8*node+1];                                        \
    w = tab_w[16*node];                                             \
    wpre = tab_w[16*node+1];                                        \
    Iw = tab_w[16*node+2];                                          \
    Iwpre = tab_w[16*node+3];                                       \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                              \
                       w2, w2pre, w, wpre, Iw, Iwpre,               \
                       n, n2);                                      \
} while(0)

#define TFT16_8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,       \
                              p8, p9, p10, p11, p12, p13, p14, p15, \
                              node, n, n2, tab_w)                   \
do {                                                                \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                            \
                                                                    \
    w2 = tab_w[2*node];                                             \
    w2pre = tab_w[2*node+1];                                        \
    w = tab_w[4*node];                                              \
    wpre = tab_w[4*node+1];                                         \
                                                                    \
    TFT4_2_NODE_LAZY_4_4(p0, p4, p8, p12,                           \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p1, p5, p9, p13,                           \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p2, p6, p10, p14,                          \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
    TFT4_2_NODE_LAZY_4_4(p3, p7, p11, p15,                          \
                         w2, w2pre, w, wpre,                        \
                         n, n2);                                    \
                                                                    \
    w2 = tab_w[8*node];                                             \
    w2pre = tab_w[8*node+1];                                        \
    w = tab_w[16*node];                                             \
    wpre = tab_w[16*node+1];                                        \
    Iw = tab_w[16*node+2];                                          \
    Iwpre = tab_w[16*node+3];                                       \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                              \
                       w2, w2pre, w, wpre, Iw, Iwpre,               \
                       n, n2);                                      \
                                                                    \
    w2 = tab_w[8*node+2];                                           \
    w2pre = tab_w[8*node+3];                                        \
    w = tab_w[16*node+4];                                           \
    wpre = tab_w[16*node+5];                                        \
    Iw = tab_w[16*node+6];                                          \
    Iwpre = tab_w[16*node+7];                                       \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                              \
                       w2, w2pre, w, wpre, Iw, Iwpre,               \
                       n, n2);                                      \
} while(0)

#define TFT16_12_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,       \
                               p8, p9, p10, p11, p12, p13, p14, p15, \
                               node, n, n2, tab_w)                   \
do {                                                                 \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                             \
                                                                     \
    w2 = tab_w[2*node];                                              \
    w2pre = tab_w[2*node+1];                                         \
    w = tab_w[4*node];                                               \
    wpre = tab_w[4*node+1];                                          \
    Iw = tab_w[4*node+2];                                            \
    Iwpre = tab_w[4*node+3];                                         \
                                                                     \
    DFT4_NODE_LAZY_4_4(p0, p4, p8, p12,                              \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
    DFT4_NODE_LAZY_4_4(p1, p5, p9, p13,                              \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
    DFT4_NODE_LAZY_4_4(p2, p6, p10, p14,                             \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
    DFT4_NODE_LAZY_4_4(p3, p7, p11, p15,                             \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
                                                                     \
    w2 = tab_w[8*node];                                              \
    w2pre = tab_w[8*node+1];                                         \
    w = tab_w[16*node];                                              \
    wpre = tab_w[16*node+1];                                         \
    Iw = tab_w[16*node+2];                                           \
    Iwpre = tab_w[16*node+3];                                        \
    DFT4_NODE_LAZY_4_4(p0, p1, p2, p3,                               \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
                                                                     \
    w2 = tab_w[8*node+2];                                            \
    w2pre = tab_w[8*node+3];                                         \
    w = tab_w[16*node+4];                                            \
    wpre = tab_w[16*node+5];                                         \
    Iw = tab_w[16*node+6];                                           \
    Iwpre = tab_w[16*node+7];                                        \
    DFT4_NODE_LAZY_4_4(p4, p5, p6, p7,                               \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
                                                                     \
    w2 = tab_w[8*node+4];                                            \
    w2pre = tab_w[8*node+5];                                         \
    w = tab_w[16*node+8];                                            \
    wpre = tab_w[16*node+9];                                         \
    Iw = tab_w[16*node+10];                                          \
    Iwpre = tab_w[16*node+11];                                       \
    DFT4_NODE_LAZY_4_4(p8, p9, p10, p11,                             \
                       w2, w2pre, w, wpre, Iw, Iwpre,                \
                       n, n2);                                       \
} while(0)

/*-------------------------*/
/* length 32, general node */
/*-------------------------*/

/** 32-point FFT, evaluation, general node, truncated at {4,8,12,16,20,24,28}
 * * same as DFT32_NODE_LAZY_4_4 but only computes first k evaluations (for k
 * in {4, 8, 12, 16, 20, 24, 28}, 7 variants available) in the first k elements
 * of p; other 32 - k elements may be modified
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
#define TFT32_4_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                    \
                              p8, p9, p10, p11, p12, p13, p14, p15,              \
                              p16, p17, p18, p19, p20, p21, p22, p23,            \
                              p24, p25, p26, p27, p28, p29, p30, p31,            \
                              node, n, n2, tab_w)                                \
do {                                                                             \
    ulong w2 = tab_w[2*node];                                                    \
    ulong w2pre = tab_w[2*node+1];                                               \
    ulong w = tab_w[4*node];                                                     \
    ulong wpre = tab_w[4*node+1];                                                \
    TFT4_2_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, n, n2);           \
    TFT4_2_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, n, n2);           \
    TFT4_2_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, n, n2);          \
    TFT4_2_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, n, n2);          \
    TFT4_2_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, n, n2);          \
    TFT4_2_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, n, n2);          \
    TFT4_2_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, n, n2);          \
    TFT4_2_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, n, n2);          \
                                                                                 \
    TFT8_4_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);  \
} while(0)

#define TFT32_8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                  \
                              p8, p9, p10, p11, p12, p13, p14, p15,            \
                              p16, p17, p18, p19, p20, p21, p22, p23,          \
                              p24, p25, p26, p27, p28, p29, p30, p31,          \
                              node, n, n2, tab_w)                              \
do {                                                                           \
    ulong w2 = tab_w[2*node];                                                  \
    ulong w2pre = tab_w[2*node+1];                                             \
    ulong w = tab_w[4*node];                                                   \
    ulong wpre = tab_w[4*node+1];                                              \
    TFT4_2_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, n, n2);         \
    TFT4_2_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, n, n2);         \
    TFT4_2_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, n, n2);        \
    TFT4_2_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, n, n2);        \
    TFT4_2_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, n, n2);        \
    TFT4_2_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, n, n2);        \
    TFT4_2_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, n, n2);        \
    TFT4_2_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, n, n2);        \
                                                                               \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);  \
} while(0)

#define TFT32_12_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                           \
                               p8, p9, p10, p11, p12, p13, p14, p15,                     \
                               p16, p17, p18, p19, p20, p21, p22, p23,                   \
                               p24, p25, p26, p27, p28, p29, p30, p31,                   \
                               node, n, n2, tab_w)                                       \
do {                                                                                     \
    ulong w2 = tab_w[2*node];                                                            \
    ulong w2pre = tab_w[2*node+1];                                                       \
    ulong w = tab_w[4*node];                                                             \
    ulong wpre = tab_w[4*node+1];                                                        \
    TFT4_2_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, n, n2);                   \
    TFT4_2_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, n, n2);                   \
    TFT4_2_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, n, n2);                  \
    TFT4_2_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, n, n2);                  \
    TFT4_2_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, n, n2);                  \
    TFT4_2_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, n, n2);                  \
    TFT4_2_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, n, n2);                  \
    TFT4_2_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, n, n2);                  \
                                                                                         \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);            \
    TFT8_4_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);  \
} while(0)

#define TFT32_16_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                         \
                               p8, p9, p10, p11, p12, p13, p14, p15,                   \
                               p16, p17, p18, p19, p20, p21, p22, p23,                 \
                               p24, p25, p26, p27, p28, p29, p30, p31,                 \
                               node, n, n2, tab_w)                                     \
do {                                                                                   \
    ulong w2 = tab_w[2*node];                                                          \
    ulong w2pre = tab_w[2*node+1];                                                     \
    ulong w = tab_w[4*node];                                                           \
    ulong wpre = tab_w[4*node+1];                                                      \
    TFT4_2_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, n, n2);                 \
    TFT4_2_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, n, n2);                 \
    TFT4_2_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, n, n2);                \
    TFT4_2_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, n, n2);                \
    TFT4_2_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, n, n2);                \
    TFT4_2_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, n, n2);                \
    TFT4_2_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, n, n2);                \
    TFT4_2_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, n, n2);                \
                                                                                       \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);          \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);  \
} while(0)

#define TFT32_20_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                             \
                               p8, p9, p10, p11, p12, p13, p14, p15,                       \
                               p16, p17, p18, p19, p20, p21, p22, p23,                     \
                               p24, p25, p26, p27, p28, p29, p30, p31,                     \
                               node, n, n2, tab_w)                                         \
do {                                                                                       \
    ulong w2 = tab_w[2*node];                                                              \
    ulong w2pre = tab_w[2*node+1];                                                         \
    ulong w = tab_w[4*node];                                                               \
    ulong wpre = tab_w[4*node+1];                                                          \
    ulong Iw = tab_w[4*node+2];                                                            \
    ulong Iwpre = tab_w[4*node+3];                                                         \
    DFT4_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);            \
    DFT4_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);            \
    DFT4_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
                                                                                           \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);              \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);      \
    TFT8_4_NODE_LAZY_4_4(p16, p17, p18, p19, p20, p21, p22, p23, 4*node+2, n, n2, tab_w);  \
} while(0)

#define TFT32_24_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                           \
                               p8, p9, p10, p11, p12, p13, p14, p15,                     \
                               p16, p17, p18, p19, p20, p21, p22, p23,                   \
                               p24, p25, p26, p27, p28, p29, p30, p31,                   \
                               node, n, n2, tab_w)                                       \
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
} while(0)

#define TFT32_28_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7,                             \
                               p8, p9, p10, p11, p12, p13, p14, p15,                       \
                               p16, p17, p18, p19, p20, p21, p22, p23,                     \
                               p24, p25, p26, p27, p28, p29, p30, p31,                     \
                               node, n, n2, tab_w)                                         \
do {                                                                                       \
    ulong w2 = tab_w[2*node];                                                              \
    ulong w2pre = tab_w[2*node+1];                                                         \
    ulong w = tab_w[4*node];                                                               \
    ulong wpre = tab_w[4*node+1];                                                          \
    ulong Iw = tab_w[4*node+2];                                                            \
    ulong Iwpre = tab_w[4*node+3];                                                         \
    DFT4_NODE_LAZY_4_4(p0, p8, p16, p24, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);            \
    DFT4_NODE_LAZY_4_4(p1, p9, p17, p25, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);            \
    DFT4_NODE_LAZY_4_4(p2, p10, p18, p26, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p3, p11, p19, p27, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p4, p12, p20, p28, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p5, p13, p21, p29, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p6, p14, p22, p30, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
    DFT4_NODE_LAZY_4_4(p7, p15, p23, p31, w2, w2pre, w, wpre, Iw, Iwpre, n, n2);           \
                                                                                           \
    DFT8_NODE_LAZY_4_4(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, n, n2, tab_w);              \
    DFT8_NODE_LAZY_4_4(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, n, n2, tab_w);      \
    DFT8_NODE_LAZY_4_4(p16, p17, p18, p19, p20, p21, p22, p23, 4*node+2, n, n2, tab_w);  \
    TFT8_4_NODE_LAZY_4_4(p24, p25, p26, p27, p28, p29, p30, p31, 4*node+3, n, n2, tab_w);  \
} while(0)


#endif  /* N_FFT_MACROS_TFT_H */

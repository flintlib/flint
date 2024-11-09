/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft.h"
#include "n_fft_macros.h"

/*-------------*/
/* 2-point DFT */
/*-------------*/

/** Cooley-Tukey butterfly, general
 * * in [0..4n) / out [0..4n) / max < 4n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [w -w]
 * * n2 is 2*n, w_pr is the precomputed data for multiplication by w mod n
 *   p_hi, p_lo, u, v are temporaries
 */
#define DFT2_LAZY44(a, b, n, n2, w, w_pr, p_hi, p_lo, u, v) \
    do {                                                          \
        u = (a);                                                  \
        if (u >= (n2))                                            \
            u -= (n2);  /* [0..2n) */                             \
        v = (b);                                                  \
        N_MULMOD_PRECOMP_LAZY(v, w, v, w_pr, n, p_hi, p_lo);      \
        (a) = u + v;                   /* [0..4n) */              \
        (b) = u + (n2) - v;         /* [0..4n) */                 \
    } while(0)

/*-------------*/
/* 4-point DFT */
/*-------------*/

/** 4-point DFT, general
 * * in [0..4n) / out [0..4n) / max < 4n
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
 */
#define DFT4_LAZY44(a, b, c, d,                                \
                   w1, w1_pr, w2, w2_pr, w3, w3_pr,            \
                   n, n2, p_hi, p_lo, tmp)                     \
do {                                                           \
    ulong u0 = (a);                                            \
    ulong u1 = (b);                                            \
    ulong u2 = (c);                                            \
    ulong u3 = (d);                                            \
    if (u0 >= n2)                                              \
        u0 -= n2;                                              \
    if (u1 >= n2)                                              \
        u1 -= n2;                                              \
                                                               \
    N_MULMOD_PRECOMP_LAZY(u2, w1, u2, w1_pr, n, p_hi, p_lo);   \
    tmp = u0;                                                  \
    u0 = u0 + u2;                    /* [0..4n) */             \
    u2 = tmp + n2 - u2;              /* [0..4n) */             \
    if (u0 >= n2)                                              \
        u0 -= n2;                    /* [0..2n) */             \
    if (u2 >= n2)                                              \
        u2 -= n2;                    /* [0..2n) */             \
                                                               \
    N_MULMOD_PRECOMP_LAZY(u3, w1, u3, w1_pr, n, p_hi, p_lo);   \
    tmp = u1;                                                  \
    u1 = u1 + u3;                    /* [0..4n) */             \
    u3 = tmp + n2 - u3;              /* [0..4n) */             \
                                                               \
    N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n, p_hi, p_lo);   \
    tmp = u0;                                                  \
    (a) = u0 + u1;                   /* [0..4n) */             \
    (b) = tmp + n2 - u1;             /* [0..4n) */             \
                                                               \
    N_MULMOD_PRECOMP_LAZY(u3, w3, u3, w3_pr, n, p_hi, p_lo);   \
    tmp = u2;                                                  \
    (c) = u2 + u3;                    /* [0..4n) */            \
    (d) = tmp + n2 - u3;              /* [0..4n) */            \
} while(0)

/*-------------*/
/* 8-point DFT */
/*-------------*/

/** 8-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 4n
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7 into its evaluations
 *       p(1), p(-1), p(I), p(-I), p(J), p(-J), p(I*J), p(-I*J)
 * i.e. the evaluations at all 8-th roots of unity J**k for 0 <= k < 8 in
 * bit-reversed order
 * * Recall [F->tab_w[2*k] for k in range(4)] == [1, I, J, IJ]
 */
#define DFT8_NODE0_LAZY14(p0, p1, p2, p3, p4, p5, p6, p7,  \
                          mod, mod2, tab_w)                \
do {                                                       \
    ulong p_hi, p_lo, tmp;                                 \
                                                           \
    BUTTERFLY_LAZY12(p0, p4, mod, tmp);                    \
    BUTTERFLY_LAZY12(p1, p5, mod, tmp);                    \
    BUTTERFLY_LAZY12(p2, p6, mod, tmp);                    \
    BUTTERFLY_LAZY12(p3, p7, mod, tmp);                    \
                                                           \
    DFT4_NODE0_LAZY24(p0, p1, p2, p3,                      \
                      tab_w[2], tab_w[3],                  \
                      mod, mod2, p_hi, p_lo);              \
    /* could use a lazy24 variant of the next macro, */    \
    /* but the gain is negligible                    */    \
    DFT4_LAZY44(p4, p5, p6, p7,                            \
                tab_w[2], tab_w[3],                        \
                tab_w[4], tab_w[5],                        \
                tab_w[6], tab_w[7],                        \
                mod, mod2, p_hi, p_lo, tmp);               \
} while(0)

/** 8-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 4n
 * * apart from these ranges, same specification as DFT8_NODE0_LAZY14
 */
#define DFT8_NODE0_LAZY24(p0, p1, p2, p3, p4, p5, p6, p7,  \
                          mod, mod2, tab_w)                \
do {                                                       \
    ulong p_hi, p_lo, tmp;                                 \
                                                           \
    BUTTERFLY_LAZY24(p0, p4, mod2, tmp);                   \
    BUTTERFLY_LAZY24(p1, p5, mod2, tmp);                   \
    BUTTERFLY_LAZY24(p2, p6, mod2, tmp);                   \
    BUTTERFLY_LAZY24(p3, p7, mod2, tmp);                   \
                                                           \
    DFT4_NODE0_LAZY24(p0, p1, p2, p3,                      \
                      tab_w[2], tab_w[3],                  \
                      mod, mod2, p_hi, p_lo);              \
    DFT4_LAZY44(p4, p5, p6, p7,                            \
                tab_w[2], tab_w[3],                        \
                tab_w[4], tab_w[5],                        \
                tab_w[6], tab_w[7],                        \
                mod, mod2, p_hi, p_lo, tmp);               \
} while(0)

/** 8-point DFT
 * * in [0..4n) / out [0..4n) / max < 4n
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7 into its evaluations
 *       p(w0), p(-w0), p(w1), p(-w1), p(w2), p(-w2), p(w3), p(-w3)
 * where w_k = F->tab_w[8*node + 2*k] for 0 <= k < 4
 * * By construction these 8 evaluation points are the 8 roots of the
 * polynomial x**8 - F->tab_w[node]
 */
#define DFT8_LAZY44(p0, p1, p2, p3, p4, p5, p6, p7,              \
                    node, mod, mod2, tab_w)                      \
do {                                                             \
    ulong p_hi, p_lo, u, v;                                      \
                                                                 \
    const ulong w = tab_w[2*(node)];                             \
    const ulong w_pr = tab_w[2*(node)+1];                        \
    DFT2_LAZY44(p0, p4, mod, mod2, w, w_pr, p_hi, p_lo, u, v);   \
    DFT2_LAZY44(p1, p5, mod, mod2, w, w_pr, p_hi, p_lo, u, v);   \
    DFT2_LAZY44(p2, p6, mod, mod2, w, w_pr, p_hi, p_lo, u, v);   \
    DFT2_LAZY44(p3, p7, mod, mod2, w, w_pr, p_hi, p_lo, u, v);   \
                                                                 \
    DFT4_LAZY44(p0, p1, p2, p3,                                  \
                tab_w[4*(node)], tab_w[4*(node)+1],              \
                tab_w[8*(node)], tab_w[8*(node)+1],              \
                tab_w[8*(node)+2], tab_w[8*(node)+3],            \
                mod, mod2, p_hi, p_lo, u);                       \
                                                                 \
    DFT4_LAZY44(p4, p5, p6, p7,                                  \
                tab_w[4*(node)+2], tab_w[4*(node)+3],            \
                tab_w[8*(node)+4], tab_w[8*(node)+5],            \
                tab_w[8*(node)+6], tab_w[8*(node)+7],            \
                mod, mod2, p_hi, p_lo, u);                       \
} while(0)

/*--------------*/
/* 16-point DFT */
/*--------------*/

/** 16-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 4n
 * * Apart from this range, same specification as dft_node0_lazy24, for depth==4
 */
#define DFT16_NODE0_LAZY14(p0, p1, p2, p3, p4, p5, p6, p7,        \
                           p8, p9, p10, p11, p12, p13, p14, p15,  \
                           mod, mod2, tab_w)                      \
do {                                                              \
    ulong p_hi, p_lo, tmp;                                        \
                                                                  \
    DFT4_NODE0_LAZY14(p0, p4, p8, p12,                            \
                      tab_w[2], tab_w[3],                         \
                      mod, mod2, p_hi, p_lo);                     \
    if (p0 >= mod2)                                               \
        p0 -= mod2;                                               \
    DFT4_NODE0_LAZY14(p1, p5, p9, p13,                            \
                      tab_w[2], tab_w[3],                         \
                      mod, mod2, p_hi, p_lo);                     \
    if (p1 >= mod2)                                               \
        p1 -= mod2;                                               \
    DFT4_NODE0_LAZY14(p2, p6, p10, p14,                           \
                      tab_w[2], tab_w[3],                         \
                      mod, mod2, p_hi, p_lo);                     \
    if (p2 >= mod2)                                               \
        p2 -= mod2;                                               \
    DFT4_NODE0_LAZY14(p3, p7, p11, p15,                           \
                      tab_w[2], tab_w[3],                         \
                      mod, mod2, p_hi, p_lo);                     \
    if (p3 >= mod2)                                               \
        p3 -= mod2;                                               \
                                                                  \
    /* next line requires < 2n,        */                         \
    /* hence the four reductions above */                         \
    DFT4_NODE0_LAZY24(p0, p1, p2, p3,                             \
                      tab_w[2], tab_w[3],                         \
                      mod, mod2, p_hi, p_lo);                     \
    DFT4_LAZY44(p4, p5, p6, p7,                                   \
                tab_w[2], tab_w[3],                               \
                tab_w[4], tab_w[5],                               \
                tab_w[6], tab_w[7],                               \
                mod, mod2, p_hi, p_lo, tmp);                      \
    DFT4_LAZY44(p8, p9, p10, p11,                                 \
                tab_w[4], tab_w[5],                               \
                tab_w[8], tab_w[9],                               \
                tab_w[10], tab_w[11],                             \
                mod, mod2, p_hi, p_lo, tmp);                      \
    DFT4_LAZY44(p12, p13, p14, p15,                               \
                tab_w[6], tab_w[7],                               \
                tab_w[12], tab_w[13],                             \
                tab_w[14], tab_w[15],                             \
                mod, mod2, p_hi, p_lo, tmp);                      \
} while(0)

/** 16-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 4n
 * * Same specification as dft_node0_lazy24, for depth==4
 */
#define DFT16_NODE0_LAZY24(p0, p1, p2, p3, p4, p5, p6, p7,       \
                           p8, p9, p10, p11, p12, p13, p14, p15, \
                           mod, mod2, tab_w)                     \
do {                                                             \
    ulong p_hi, p_lo, tmp;                                       \
                                                                 \
    DFT4_NODE0_LAZY24(p0, p4, p8, p12,                           \
                      tab_w[2], tab_w[3],                        \
                      mod, mod2, p_hi, p_lo);                    \
    if (p0 >= mod2)                                              \
        p0 -= mod2;                                              \
    DFT4_NODE0_LAZY24(p1, p5, p9, p13,                           \
                      tab_w[2], tab_w[3],                        \
                      mod, mod2, p_hi, p_lo);                    \
    if (p1 >= mod2)                                              \
        p1 -= mod2;                                              \
    DFT4_NODE0_LAZY24(p2, p6, p10, p14,                          \
                      tab_w[2], tab_w[3],                        \
                      mod, mod2, p_hi, p_lo);                    \
    if (p2 >= mod2)                                              \
        p2 -= mod2;                                              \
    DFT4_NODE0_LAZY24(p3, p7, p11, p15,                          \
                      tab_w[2], tab_w[3],                        \
                      mod, mod2, p_hi, p_lo);                    \
    if (p3 >= mod2)                                              \
        p3 -= mod2;                                              \
                                                                 \
    /* next line requires < 2n,        */                        \
    /* hence the four reductions above */                        \
    DFT4_NODE0_LAZY24(p0, p1, p2, p3,                            \
                      tab_w[2], tab_w[3],                        \
                      mod, mod2, p_hi, p_lo);                    \
    DFT4_LAZY44(p4, p5, p6, p7,                                  \
                tab_w[2], tab_w[3],                              \
                tab_w[4], tab_w[5],                              \
                tab_w[6], tab_w[7],                              \
                mod, mod2, p_hi, p_lo, tmp);                     \
    DFT4_LAZY44(p8, p9, p10, p11,                                \
                tab_w[4], tab_w[5],                              \
                tab_w[8], tab_w[9],                              \
                tab_w[10], tab_w[11],                            \
                mod, mod2, p_hi, p_lo, tmp);                     \
    DFT4_LAZY44(p12, p13, p14, p15,                              \
                tab_w[6], tab_w[7],                              \
                tab_w[12], tab_w[13],                            \
                tab_w[14], tab_w[15],                            \
                mod, mod2, p_hi, p_lo, tmp);                     \
} while(0)

/** 16-point DFT
 * * in [0..4n) / out [0..4n) / max < 4n
 * * Same specification as dft_lazy44, for depth==4
 */
#define DFT16_LAZY44(p0, p1, p2, p3, p4, p5, p6, p7,       \
                     p8, p9, p10, p11, p12, p13, p14, p15, \
                     node, mod, mod2, tab_w)               \
do {                                                       \
    ulong p_hi, p_lo, tmp;                                 \
    ulong w2, w2pre, w, wpre, Iw, Iwpre;                   \
                                                           \
    w2 = tab_w[2*node];                                    \
    w2pre = tab_w[2*node+1];                               \
    w = tab_w[4*node];                                     \
    wpre = tab_w[4*node+1];                                \
    Iw = tab_w[4*node+2];                                  \
    Iwpre = tab_w[4*node+3];                               \
                                                           \
    DFT4_LAZY44(p0, p4, p8, p12,                           \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
    DFT4_LAZY44(p1, p5, p9, p13,                           \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
    DFT4_LAZY44(p2, p6, p10, p14,                          \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
    DFT4_LAZY44(p3, p7, p11, p15,                          \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
                                                           \
    w2 = tab_w[8*node];                                    \
    w2pre = tab_w[8*node+1];                               \
    w = tab_w[16*node];                                    \
    wpre = tab_w[16*node+1];                               \
    Iw = tab_w[16*node+2];                                 \
    Iwpre = tab_w[16*node+3];                              \
    DFT4_LAZY44(p0, p1, p2, p3,                            \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
                                                           \
    w2 = tab_w[8*node+2];                                  \
    w2pre = tab_w[8*node+3];                               \
    w = tab_w[16*node+4];                                  \
    wpre = tab_w[16*node+5];                               \
    Iw = tab_w[16*node+6];                                 \
    Iwpre = tab_w[16*node+7];                              \
    DFT4_LAZY44(p4, p5, p6, p7,                            \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
                                                           \
    w2 = tab_w[8*node+4];                                  \
    w2pre = tab_w[8*node+5];                               \
    w = tab_w[16*node+8];                                  \
    wpre = tab_w[16*node+9];                               \
    Iw = tab_w[16*node+10];                                \
    Iwpre = tab_w[16*node+11];                             \
    DFT4_LAZY44(p8, p9, p10, p11,                          \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
                                                           \
    w2 = tab_w[8*node+6];                                  \
    w2pre = tab_w[8*node+7];                               \
    w = tab_w[16*node+12];                                 \
    wpre = tab_w[16*node+13];                              \
    Iw = tab_w[16*node+14];                                \
    Iwpre = tab_w[16*node+15];                             \
    DFT4_LAZY44(p12, p13, p14, p15,                        \
                w2, w2pre, w, wpre, Iw, Iwpre,             \
                mod, mod2, p_hi, p_lo, tmp);               \
} while(0)

/*--------------*/
/* 32-point DFT */
/*--------------*/

/** 32-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 4n
 * * Apart from this range, same specification as dft_node0_lazy24, for depth==5
 */
#define DFT32_NODE0_LAZY14(p0, p1, p2, p3, p4, p5, p6, p7,                    \
                           p8, p9, p10, p11, p12, p13, p14, p15,              \
                           p16, p17, p18, p19, p20, p21, p22, p23,            \
                           p24, p25, p26, p27, p28, p29, p30, p31,            \
                           mod, mod2, tab_w)                                  \
do {                                                                          \
    ulong p_hi, p_lo;                                                         \
                                                                              \
    DFT4_NODE0_LAZY14(p0, p8, p16, p24,                                       \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p0 >= mod2)                                                           \
        p0 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p1, p9, p17, p25,                                       \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p1 >= mod2)                                                           \
        p1 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p2, p10, p18, p26,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p2 >= mod2)                                                           \
        p2 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p3, p11, p19, p27,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p3 >= mod2)                                                           \
        p3 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p4, p12, p20, p28,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p4 >= mod2)                                                           \
        p4 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p5, p13, p21, p29,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p5 >= mod2)                                                           \
        p5 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p6, p14, p22, p30,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p6 >= mod2)                                                           \
        p6 -= mod2;                                                           \
    DFT4_NODE0_LAZY14(p7, p15, p23, p31,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p7 >= mod2)                                                           \
        p7 -= mod2;                                                           \
                                                                              \
    /* next line requires < 2n, hence the 8 reductions above */               \
    DFT8_NODE0_LAZY24(p0, p1, p2, p3, p4, p5, p6, p7, mod, mod2, tab_w);      \
    DFT8_LAZY44(p8, p9, p10, p11, p12, p13, p14, p15, 1, mod, mod2, tab_w);   \
    DFT8_LAZY44(p16, p17, p18, p19, p20, p21, p22, p23, 2, mod, mod2, tab_w); \
    DFT8_LAZY44(p24, p25, p26, p27, p28, p29, p30, p31, 3, mod, mod2, tab_w); \
} while(0)

/** 32-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 4n
 * * Same specification as dft_node0_lazy24, for depth==5
 */
#define DFT32_NODE0_LAZY24(p0, p1, p2, p3, p4, p5, p6, p7,                    \
                           p8, p9, p10, p11, p12, p13, p14, p15,              \
                           p16, p17, p18, p19, p20, p21, p22, p23,            \
                           p24, p25, p26, p27, p28, p29, p30, p31,            \
                           mod, mod2, tab_w)                                  \
do {                                                                          \
    ulong p_hi, p_lo;                                                         \
                                                                              \
    DFT4_NODE0_LAZY24(p0, p8, p16, p24,                                       \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p0 >= mod2)                                                           \
        p0 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p1, p9, p17, p25,                                       \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p1 >= mod2)                                                           \
        p1 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p2, p10, p18, p26,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p2 >= mod2)                                                           \
        p2 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p3, p11, p19, p27,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p3 >= mod2)                                                           \
        p3 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p4, p12, p20, p28,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p4 >= mod2)                                                           \
        p4 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p5, p13, p21, p29,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p5 >= mod2)                                                           \
        p5 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p6, p14, p22, p30,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p6 >= mod2)                                                           \
        p6 -= mod2;                                                           \
    DFT4_NODE0_LAZY24(p7, p15, p23, p31,                                      \
                      tab_w[2], tab_w[3],                                     \
                      mod, mod2, p_hi, p_lo);                                 \
    if (p7 >= mod2)                                                           \
        p7 -= mod2;                                                           \
                                                                              \
    /* next line requires < 2n, hence the 8 reductions above */               \
    DFT8_NODE0_LAZY24(p0, p1, p2, p3, p4, p5, p6, p7, mod, mod2, tab_w);      \
    DFT8_LAZY44(p8, p9, p10, p11, p12, p13, p14, p15, 1, mod, mod2, tab_w);   \
    DFT8_LAZY44(p16, p17, p18, p19, p20, p21, p22, p23, 2, mod, mod2, tab_w); \
    DFT8_LAZY44(p24, p25, p26, p27, p28, p29, p30, p31, 3, mod, mod2, tab_w); \
} while(0)

/** 32-point DFT
 * * in [0..4n) / out [0..4n) / max < 4n
 * * Same specification as dft_lazy44, for depth==5
 */
#define DFT32_LAZY44(p0, p1, p2, p3, p4, p5, p6, p7,                                           \
                     p8, p9, p10, p11, p12, p13, p14, p15,                                     \
                     p16, p17, p18, p19, p20, p21, p22, p23,                                   \
                     p24, p25, p26, p27, p28, p29, p30, p31,                                   \
                     node, mod, mod2, tab_w)                                                   \
do {                                                                                           \
    ulong p_hi, p_lo, tmp;                                                                     \
                                                                                               \
    ulong w2 = tab_w[2*node];                                                                  \
    ulong w2pre = tab_w[2*node+1];                                                             \
    ulong w = tab_w[4*node];                                                                   \
    ulong wpre = tab_w[4*node+1];                                                              \
    ulong Iw = tab_w[4*node+2];                                                                \
    ulong Iwpre = tab_w[4*node+3];                                                             \
    DFT4_LAZY44(p0, p8, p16, p24, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp);  \
    DFT4_LAZY44(p1, p9, p17, p25, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp);  \
    DFT4_LAZY44(p2, p10, p18, p26, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
    DFT4_LAZY44(p3, p11, p19, p27, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
    DFT4_LAZY44(p4, p12, p20, p28, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
    DFT4_LAZY44(p5, p13, p21, p29, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
    DFT4_LAZY44(p6, p14, p22, p30, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
    DFT4_LAZY44(p7, p15, p23, p31, w2, w2pre, w, wpre, Iw, Iwpre, mod, mod2, p_hi, p_lo, tmp); \
                                                                                               \
    /* next line requires < 2n, hence the four reductions above */                             \
    DFT8_LAZY44(p0, p1, p2, p3, p4, p5, p6, p7, 4*node, mod, mod2, tab_w);                     \
    DFT8_LAZY44(p8, p9, p10, p11, p12, p13, p14, p15, 4*node+1, mod, mod2, tab_w);             \
    DFT8_LAZY44(p16, p17, p18, p19, p20, p21, p22, p23, 4*node+2, mod, mod2, tab_w);           \
    DFT8_LAZY44(p24, p25, p26, p27, p28, p29, p30, p31, 4*node+3, mod, mod2, tab_w);           \
} while(0)

/*-------------*/
/* general DFT */
/*-------------*/

/** 2**depth-point DFT
 * * in [0..4n) / out [0..4n) / max < 4n
 * * In-place transform p of length len == 2**depth into
 * the concatenation of
 *       [sum(p[i] * w_k**i for i in range(len), sum(p[i] * (-w_k)**i for i in range(len)]
 * for k in range(len),
 * where w_k = F->tab_w[2**depth * node + 2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the roots of the
 * polynomial x**len - F->tab_w[node]
 * * Requirement (not checked):
 *        3 <= depth
 *        (node+1) * 2**depth <= 2**F.depth (length of F->tab_w)
 */
void dft_lazy44(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 3)
    {
        DFT8_LAZY44(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        DFT16_LAZY44(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                     p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                     node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_LAZY44(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                     p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                     p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                     p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                     node, F->mod, F->mod2, F->tab_w);
    }
    else
    {
        const ulong len = UWORD(1) << depth;

        // 4-point butterflies
        // in: [0..4n), out: [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p+len/4;
        const nn_ptr p2 = p+2*len/4;
        const nn_ptr p3 = p+3*len/4;
        const ulong w2 = F->tab_w[2*node];
        const ulong w2pre = F->tab_w[2*node+1];
        const ulong w = F->tab_w[4*node];
        const ulong wpre = F->tab_w[4*node+1];
        const ulong Iw = F->tab_w[4*node+2];
        const ulong Iwpre = F->tab_w[4*node+3];
        ulong p_hi, p_lo, tmp;

        for (ulong k = 0; k < len/4; k+=4)
        {
            DFT4_LAZY44(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY44(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY44(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_LAZY44(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
        }

        // 4 recursive calls with depth-2
        dft_lazy44(p0, depth-2, 4*node, F);
        dft_lazy44(p1, depth-2, 4*node+1, F);
        dft_lazy44(p2, depth-2, 4*node+2, F);
        dft_lazy44(p3, depth-2, 4*node+3, F);
    }
}

/** 2**depth-point DFT
 * * in [0..2n) / out [0..4n) / max < 4n
 * * In-place transform p of length len == 2**depth into
 * the concatenation of
 *       [sum(p[i] * w_k**i for i in range(len), sum(p[i] * (-w_k)**i for i in range(len)]
 * for k in range(len),
 * where w_k = F->tab_w[2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the roots of the polynomial
 * x**len - 1, precisely they are all powers of the chosen len-th primitive
 * root of unity with exponents listed in bit reversed order
 * * Requirements (not checked): 3 <= depth <= F.depth
 */
void dft_node0_lazy24(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 3)
    {
        DFT8_NODE0_LAZY24(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        DFT16_NODE0_LAZY24(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                           p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                           F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_NODE0_LAZY24(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                           p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                           p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                           p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                           F->mod, F->mod2, F->tab_w);
    }
    else
    {
        const ulong len = UWORD(1) << depth;

        // 4-point butterflies
        // input p0,p1,p2,p3 in [0..2n) x [0..2n) x [0..2n) x [0..2n)
        // output p0,p1,p2,p3 in [0..2n) x [0..4n) x [0..4n) x [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/4;
        const nn_ptr p2 = p + 2*len/4;
        const nn_ptr p3 = p + 3*len/4;
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k++)
        {
            DFT4_NODE0_LAZY24(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            if (p0[k] >= F->mod2)
                p0[k] -= F->mod2;
        }

        // 4 recursive calls with depth-2
        dft_node0_lazy24(p0, depth-2, F);
        dft_lazy44(p1, depth-2, 1, F);
        dft_lazy44(p2, depth-2, 2, F);
        dft_lazy44(p3, depth-2, 3, F);
    }
}

/** 2**depth-point DFT
 * * in [0..n) / out [0..4n) / max < 4n
 * * In-place transform p of length len == 2**depth into
 * the concatenation of
 *       [sum(p[i] * w_k**i for i in range(len), sum(p[i] * (-w_k)**i for i in range(len)]
 * for k in range(len),
 * where w_k = F->tab_w[2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the roots of the polynomial
 * x**len - 1, precisely they are all powers of the chosen len-th primitive
 * root of unity with exponents listed in bit reversed order
 * * Requirements (not checked): depth <= F.depth
 */
void dft_node0_lazy14(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 4)
    {
        DFT16_NODE0_LAZY14(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                           p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                           F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_NODE0_LAZY14(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                           p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                           p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                           p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                           F->mod, F->mod2, F->tab_w);
    }
    else if (depth > 5)
    {
        const ulong len = UWORD(1) << depth;

        // 4-point butterflies
        // input p0,p1,p2,p3 in [0..n) x [0..n) x [0..n) x [0..n)
        // output p0,p1,p2,p3 in [0..2n) x [0..4n) x [0..4n) x [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/4;
        const nn_ptr p2 = p + 2*len/4;
        const nn_ptr p3 = p + 3*len/4;
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k++)
        {
            DFT4_NODE0_LAZY14(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            if (p0[k] >= F->mod2)
                p0[k] -= F->mod2;
        }

        // 4 recursive calls with depth-2
        dft_node0_lazy24(p0, depth-2, F);
        dft_lazy44(p1, depth-2, 1, F);
        dft_lazy44(p2, depth-2, 2, F);
        dft_lazy44(p3, depth-2, 3, F);
    }
    else if (depth == 3)
    {
        DFT8_NODE0_LAZY14(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        DFT4_NODE0_LAZY14(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 1)
    {
        ulong tmp;
        BUTTERFLY_LAZY12(p[0], p[1], F->mod, tmp);
    }
}


/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "n_fft.h"

/*---------*/
/* helpers */
/*---------*/

/** Shoup's modular multiplication with precomputation, lazy
 * (does not perform the excess correction step)
 *  --> computes either r or r+n and store it is res, where r = (a*b) % n
 *  --> a_pr is the precomputation for n, p_hi and p_lo are temporaries
 *  --> requires nbits(n) < FLINT_BITS
 */
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n, p_hi, p_lo) \
    do {                                                      \
        umul_ppmm(p_hi, p_lo, (a_pr), (b));                   \
        res = (a) * (b) - p_hi * (n);                         \
    } while(0)

/*-------------*/
/* 2-point DFT */
/*-------------*/

/** Cooley-Tukey butterfly, node 0
 * * in [0..n) x [0..n) / out [0..2n) x [0..2n) / max < 2n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n is the modulus, tmp is a temporary
 */
#define DFT2_NODE0_LAZY12(a, b, n, tmp) \
    do {                                \
        tmp = (b);                      \
        (b) = (a) + (n) - tmp;          \
        (a) = (a) + tmp;                \
    } while(0)

/** Cooley-Tukey butterfly, node 0
 * * in [0..2n) x [0..2n) / out [0..2n) x [0..4n) / max < 4n
 * * In-place transform
 *                            [1  1]
 *           [a  b] <- [a  b] [1 -1]
 * * n2 is 2*n, tmp is a temporary
 */
#define DFT2_NODE0_LAZY24(a, b, n2, tmp) \
    do {                               \
        tmp = (b);                     \
        (b) = (a) + (n2) - tmp;        \
        (a) = (a) + tmp;               \
        if ((a) >= (n2))               \
            (a) -= (n2);               \
    } while(0)

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

/** 4-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 4n
 * * In-place transform
 *                              [1  1  1  1]
 *                              [1 -1  I -I]
 * [a  b  c  d] <- [a  b  c  d] [1  1 -1 -1]
 *                              [1 -1 -I  I]
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

/** 4-point DFT, node 0
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

/** 4-point DFT, general
 * * in [0..4n) / out [0..4n) / max < 8n
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
 * typically w2**2 == w1 and w3 == I*w2 (hence w3**2 == -w1) so that this
 * really is the subproduct tree built from the four roots
 *           w2, -w2, I*w2, -I*w2   of x**4 - w1
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
    u1 = u1 + u3;                    /* [0..8n) */             \
    u3 = tmp + n2 - u3;              /* [0..8n) */             \
                                                               \
    N_MULMOD_PRECOMP_LAZY(u1, w2, u1, w2_pr, n, p_hi, p_lo);   \
    tmp = u0;                                                  \
    (a) = u0 + u1;                    /* [0..4n) */            \
    (b) = tmp + n2 - u1;              /* [0..4n) */            \
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
 * * in [0..n) / out [0..4n) / max < 8n
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7 into its evaluations 
 *       p(1), p(-1), p(I), p(-I), p(J), p(-J), p(I*J), p(-I*J)
 * i.e. the evaluations at all 8-th roots of unity J**k for 0 <= k < 8 in
 * bit-reversed order
 * * Recall [F->tab_w[2*k] for k in range(4)] == [1, I, J, IJ]
 */
FLINT_FORCE_INLINE void dft8_node0_lazy14(ulong * p0, ulong * p1, ulong * p2, ulong * p3,
                                          ulong * p4, ulong * p5, ulong * p6, ulong * p7,
                                          n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT2_NODE0_LAZY12(*p0, *p4, F->mod, tmp);
    DFT2_NODE0_LAZY12(*p1, *p5, F->mod, tmp);
    DFT2_NODE0_LAZY12(*p2, *p6, F->mod, tmp);
    DFT2_NODE0_LAZY12(*p3, *p7, F->mod, tmp);

    DFT4_NODE0_LAZY24(*p0, *p1, *p2, *p3, F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    // could use a lazy24 variant of the next macro, but the gain is negligible
    DFT4_LAZY44(*p4, *p5, *p6, *p7,
                F->tab_w[2], F->tab_w[3],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[6], F->tab_w[7],
                F->mod, F->mod2, p_hi, p_lo, tmp);
}

/** 8-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 8n
 * * apart from these ranges, same specification as dft8_node0_lazy14
 */
FLINT_FORCE_INLINE void dft8_node0_lazy24(ulong * p0, ulong * p1, ulong * p2, ulong * p3,
                                          ulong * p4, ulong * p5, ulong * p6, ulong * p7,
                                          n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT2_NODE0_LAZY24(*p0, *p4, F->mod2, tmp);
    DFT2_NODE0_LAZY24(*p1, *p5, F->mod2, tmp);
    DFT2_NODE0_LAZY24(*p2, *p6, F->mod2, tmp);
    DFT2_NODE0_LAZY24(*p3, *p7, F->mod2, tmp);

    DFT4_NODE0_LAZY24(*p0, *p1, *p2, *p3, F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY44(*p4, *p5, *p6, *p7,
                F->tab_w[2], F->tab_w[3],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[6], F->tab_w[7],
                F->mod, F->mod2, p_hi, p_lo, tmp);
}

/** 8-point DFT
 * * in [0..4n) / out [0..4n) / max < 8n
 * * In-place transform p = [p0,p1,p2,p3,p4,p5,p6,p7], seen as a polynomial
 * p(x) = p0 + p1*x + ... + p7*x**7 into its evaluations 
 *       p(w0), p(-w0), p(w1), p(-w1), p(w2), p(-w2), p(w3), p(-w3)
 * where w_k = F->tab_w[8*node + 2*k] for 0 <= k < 4
 * * By construction these 8 evaluation points are the 8 roots of the
 * polynomial x**8 - F->tab_w[node]
 */
FLINT_FORCE_INLINE void dft8_lazy44(ulong * p0, ulong * p1, ulong * p2, ulong * p3,
                                    ulong * p4, ulong * p5, ulong * p6, ulong * p7,
                                    ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, u, v;

    const ulong w = F->tab_w[2*node];
    const ulong w_pr = F->tab_w[2*node+1];
    DFT2_LAZY44(*p0, *p4, F->mod, F->mod2, w, w_pr, p_hi, p_lo, u, v);
    DFT2_LAZY44(*p1, *p5, F->mod, F->mod2, w, w_pr, p_hi, p_lo, u, v);
    DFT2_LAZY44(*p2, *p6, F->mod, F->mod2, w, w_pr, p_hi, p_lo, u, v);
    DFT2_LAZY44(*p3, *p7, F->mod, F->mod2, w, w_pr, p_hi, p_lo, u, v);

    DFT4_LAZY44(*p0, *p1, *p2, *p3,
                F->tab_w[4*node], F->tab_w[4*node+1],
                F->tab_w[8*node], F->tab_w[8*node+1],
                F->tab_w[8*node+2], F->tab_w[8*node+3],
                F->mod, F->mod2, p_hi, p_lo, u);

    DFT4_LAZY44(*p4, *p5, *p6, *p7,
                F->tab_w[4*node+2], F->tab_w[4*node+3],
                F->tab_w[8*node+4], F->tab_w[8*node+5],
                F->tab_w[8*node+6], F->tab_w[8*node+7],
                F->mod, F->mod2, p_hi, p_lo, u);
}

/*--------------*/
/* 16-point DFT */
/*--------------*/

/** 16-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 8n
 * * Apart from this range, same specification as dft_node0_lazy24, for depth==4
 */
FLINT_FORCE_INLINE void dft16_node0_lazy14(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT4_NODE0_LAZY14(p[0], p[4], p[ 8], p[12], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_NODE0_LAZY14(p[1], p[5], p[ 9], p[13], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_NODE0_LAZY14(p[2], p[6], p[10], p[14], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_NODE0_LAZY14(p[3], p[7], p[11], p[15], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;

    // next line requires < 2n, hence the four reductions above
    DFT4_NODE0_LAZY24(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY44(p[4], p[5], p[6], p[7],
                F->tab_w[2], F->tab_w[3],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[6], F->tab_w[7],
                F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[8], p[9], p[10], p[11],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[8], F->tab_w[9],
                F->tab_w[10], F->tab_w[11],
                F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[12], p[13], p[14], p[15],
                F->tab_w[6], F->tab_w[7],
                F->tab_w[12], F->tab_w[13],
                F->tab_w[14], F->tab_w[15],
                F->mod, F->mod2, p_hi, p_lo, tmp);
}

/** 16-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 8n
 * * Same specification as dft_node0_lazy24, for depth==4
 */
FLINT_FORCE_INLINE void dft16_node0_lazy24(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    DFT4_NODE0_LAZY24(p[0], p[4], p[ 8], p[12], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_NODE0_LAZY24(p[1], p[5], p[ 9], p[13], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_NODE0_LAZY24(p[2], p[6], p[10], p[14], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_NODE0_LAZY24(p[3], p[7], p[11], p[15], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;

    // next line requires < 2n, hence the four reductions above
    DFT4_NODE0_LAZY24(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    DFT4_LAZY44(p[4], p[5], p[6], p[7],
                F->tab_w[2], F->tab_w[3],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[6], F->tab_w[7],
                F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[8], p[9], p[10], p[11],
                F->tab_w[4], F->tab_w[5],
                F->tab_w[8], F->tab_w[9],
                F->tab_w[10], F->tab_w[11],
                F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[12], p[13], p[14], p[15],
                F->tab_w[6], F->tab_w[7],
                F->tab_w[12], F->tab_w[13],
                F->tab_w[14], F->tab_w[15],
                F->mod, F->mod2, p_hi, p_lo, tmp);
}

/** 16-point DFT
 * * in [0..4n) / out [0..4n) / max < 8n
 * * Same specification as dft_lazy44, for depth==4
 */
FLINT_FORCE_INLINE void dft16_lazy44(nn_ptr p, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;
    ulong w2, w2pre, w, wpre, Iw, Iwpre;

    w2 = F->tab_w[2*node];
    w2pre = F->tab_w[2*node+1];
    w = F->tab_w[4*node];
    wpre = F->tab_w[4*node+1];
    Iw = F->tab_w[4*node+2];
    Iwpre = F->tab_w[4*node+3];

    DFT4_LAZY44(p[0], p[4], p[ 8], p[12], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[1], p[5], p[ 9], p[13], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[2], p[6], p[10], p[14], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[3], p[7], p[11], p[15], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node];
    w2pre = F->tab_w[8*node+1];
    w = F->tab_w[16*node];
    wpre = F->tab_w[16*node+1];
    Iw = F->tab_w[16*node+2];
    Iwpre = F->tab_w[16*node+3];
    DFT4_LAZY44(p[0], p[1], p[2], p[3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+2];
    w2pre = F->tab_w[8*node+3];
    w = F->tab_w[16*node+4];
    wpre = F->tab_w[16*node+5];
    Iw = F->tab_w[16*node+6];
    Iwpre = F->tab_w[16*node+7];
    DFT4_LAZY44(p[4], p[5], p[6], p[7], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+4];
    w2pre = F->tab_w[8*node+5];
    w = F->tab_w[16*node+8];
    wpre = F->tab_w[16*node+9];
    Iw = F->tab_w[16*node+10];
    Iwpre = F->tab_w[16*node+11];
    DFT4_LAZY44(p[8], p[9], p[10], p[11], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    w2 = F->tab_w[8*node+6];
    w2pre = F->tab_w[8*node+7];
    w = F->tab_w[16*node+12];
    wpre = F->tab_w[16*node+13];
    Iw = F->tab_w[16*node+14];
    Iwpre = F->tab_w[16*node+15];
    DFT4_LAZY44(p[12], p[13], p[14], p[15], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
}

/*--------------*/
/* 32-point DFT */
/*--------------*/

/** 32-point DFT, node 0
 * * in [0..n) / out [0..4n) / max < 8n
 * * Apart from this range, same specification as dft_node0_lazy24, for depth==5
 */
FLINT_FORCE_INLINE void dft32_node0_lazy14(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo;

    DFT4_NODE0_LAZY14(p[0], p[8 ], p[16], p[24],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_NODE0_LAZY14(p[1], p[9 ], p[17], p[25],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_NODE0_LAZY14(p[2], p[10], p[18], p[26],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_NODE0_LAZY14(p[3], p[11], p[19], p[27],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;
    DFT4_NODE0_LAZY14(p[4], p[12], p[20], p[28],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[4] >= F->mod2)
        p[4] -= F->mod2;
    DFT4_NODE0_LAZY14(p[5], p[13], p[21], p[29],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[5] >= F->mod2)
        p[5] -= F->mod2;
    DFT4_NODE0_LAZY14(p[6], p[14], p[22], p[30],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[6] >= F->mod2)
        p[6] -= F->mod2;
    DFT4_NODE0_LAZY14(p[7], p[15], p[23], p[31],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[7] >= F->mod2)
        p[7] -= F->mod2;

    // next line requires < 2n, hence the 8 reductions above
    dft8_node0_lazy24(p+ 0, p+ 1, p+ 2, p+ 3, p+ 4, p+ 5, p+ 6, p+ 7,    F);
    dft8_lazy44(      p+ 8, p+ 9, p+10, p+11, p+12, p+13, p+14, p+15, 1, F);
    dft8_lazy44(      p+16, p+17, p+18, p+19, p+20, p+21, p+22, p+23, 2, F);
    dft8_lazy44(      p+24, p+25, p+26, p+27, p+28, p+29, p+30, p+31, 3, F);
}

/** 32-point DFT, node 0
 * * in [0..2n) / out [0..4n) / max < 8n
 * * Same specification as dft_node0_lazy24, for depth==5
 */
FLINT_FORCE_INLINE void dft32_node0_lazy24(nn_ptr p, n_fft_ctx_t F)
{
    ulong p_hi, p_lo;

    DFT4_NODE0_LAZY24(p[0], p[8 ], p[16], p[24],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[0] >= F->mod2)
        p[0] -= F->mod2;
    DFT4_NODE0_LAZY24(p[1], p[9 ], p[17], p[25],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[1] >= F->mod2)
        p[1] -= F->mod2;
    DFT4_NODE0_LAZY24(p[2], p[10], p[18], p[26],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[2] >= F->mod2)
        p[2] -= F->mod2;
    DFT4_NODE0_LAZY24(p[3], p[11], p[19], p[27],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[3] >= F->mod2)
        p[3] -= F->mod2;
    DFT4_NODE0_LAZY24(p[4], p[12], p[20], p[28],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[4] >= F->mod2)
        p[4] -= F->mod2;
    DFT4_NODE0_LAZY24(p[5], p[13], p[21], p[29],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[5] >= F->mod2)
        p[5] -= F->mod2;
    DFT4_NODE0_LAZY24(p[6], p[14], p[22], p[30],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[6] >= F->mod2)
        p[6] -= F->mod2;
    DFT4_NODE0_LAZY24(p[7], p[15], p[23], p[31],
                      F->tab_w[2], F->tab_w[3],
                      F->mod, F->mod2, p_hi, p_lo);
    if (p[7] >= F->mod2)
        p[7] -= F->mod2;

    // next line requires < 2n, hence the 8 reductions above
    dft8_node0_lazy24(p+ 0, p+ 1, p+ 2, p+ 3, p+ 4, p+ 5, p+ 6, p+ 7,    F);
    dft8_lazy44(      p+ 8, p+ 9, p+10, p+11, p+12, p+13, p+14, p+15, 1, F);
    dft8_lazy44(      p+16, p+17, p+18, p+19, p+20, p+21, p+22, p+23, 2, F);
    dft8_lazy44(      p+24, p+25, p+26, p+27, p+28, p+29, p+30, p+31, 3, F);
}

/** 32-point DFT
 * * in [0..4n) / out [0..4n) / max < 8n
 * * Same specification as dft_lazy44, for depth==5
 */
FLINT_FORCE_INLINE void dft32_lazy44(nn_ptr p, ulong node, n_fft_ctx_t F)
{
    ulong p_hi, p_lo, tmp;

    ulong w2 = F->tab_w[2*node];
    ulong w2pre = F->tab_w[2*node+1];
    ulong w = F->tab_w[4*node];
    ulong wpre = F->tab_w[4*node+1];
    ulong Iw = F->tab_w[4*node+2];
    ulong Iwpre = F->tab_w[4*node+3];
    DFT4_LAZY44(p[0], p[ 8], p[16], p[24], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[1], p[ 9], p[17], p[25], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[2], p[10], p[18], p[26], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[3], p[11], p[19], p[27], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[4], p[12], p[20], p[28], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[5], p[13], p[21], p[29], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[6], p[14], p[22], p[30], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
    DFT4_LAZY44(p[7], p[15], p[23], p[31], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);

    // next line requires < 2n, hence the four reductions above
    dft8_lazy44(p+ 0, p+ 1, p+ 2, p+ 3, p+ 4, p+ 5, p+ 6, p+ 7, 4*node, F);
    dft8_lazy44(p+ 8, p+ 9, p+10, p+11, p+12, p+13, p+14, p+15, 4*node+1, F);
    dft8_lazy44(p+16, p+17, p+18, p+19, p+20, p+21, p+22, p+23, 4*node+2, F);
    dft8_lazy44(p+24, p+25, p+26, p+27, p+28, p+29, p+30, p+31, 4*node+3, F);
}

/*-------------*/
/* general DFT */
/*-------------*/

/** 2**depth-point DFT
 * * in [0..4n) / out [0..4n) / max < 8n
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
void dft_lazy44(nn_ptr p, ulong depth, ulong node, n_fft_ctx_t F)
{
    if (depth == 3)
        dft8_lazy44(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, node, F);
    else if (depth == 4)
        dft16_lazy44(p, node, F);
    else if (depth == 5)
        dft32_lazy44(p, node, F);
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
 * * in [0..2n) / out [0..4n) / max < 8n
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
void dft_node0_lazy24(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth == 3)
        dft8_node0_lazy24(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, F);
    else if (depth == 4)
        dft16_node0_lazy24(p, F);
    else if (depth == 5)
        dft32_node0_lazy24(p, F);
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
 * * in [0..n) / out [0..4n) / max < 8n
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
void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth == 0)
        return;

    if (depth == 1)
    {
        ulong tmp;
        DFT2_NODE0_LAZY12(p[0], p[1], F->mod, tmp);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        DFT4_NODE0_LAZY14(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 3)
        dft8_node0_lazy14(p+0, p+1, p+2, p+3, p+4, p+5, p+6, p+7, F);
    else if (depth == 4)
        dft16_node0_lazy14(p, F);
    else if (depth == 5)
        dft32_node0_lazy14(p, F);
    else
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
}

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

#include "impl_macros_dft.h"
#include "n_fft/impl.h"

/*---------------------------------------*/
/* "c-circulant" division with remainder */
/*---------------------------------------*/

/* TODO bench different variants and choose fastest */

/* division by x**d - c, lazy_4_4 with precomputation */
/* in [0, 4*n) | out [0, 4*n) | max < 4n */
FLINT_FORCE_INLINE
void _nmod_poly_divrem_circulant_lazy_4_4_v0(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2)
{
    /* assumes len >= d */
    slong i;
    ulong j, r, val, p_hi, p_lo;

    r = len % d;
    i = len - r - d;  /* multiple of d, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        /* p[i+j] = p[i+j] + c * p[i+d+j] */
        val = p[d+i+j];
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        val = c * val - p_hi * n;  /* [0, 2n) */
        if (p[i+j] >= n2)
            p[i+j] -= n2;          /* [0, 2n) */
        p[i+j] = val + p[i+j];     /* [0, 4n) */
    }

    i -= d;
    while (i >= 0)
    {
        for (j = 0; j < d; j++)
        {
            /* p[i+j] = p[i+j] + c * p[i+d+j] */
            val = p[d+i+j];
            umul_ppmm(p_hi, p_lo, c_precomp, val);
            val = c * val - p_hi * n;  /* [0, 2n) */
            if (p[i+j] >= n2)
                p[i+j] -= n2;          /* [0, 2n) */
            p[i+j] = val + p[i+j];     /* [0, 4n) */
        }
        i -= d;
    }
}

/* assumes len > 0 and d > 0 */
FLINT_FORCE_INLINE
void _nmod_poly_divrem_circulant_lazy_4_4(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2)
{
    ulong i, val0, val1, p_hi, p_lo;

    for (i = len - 1; i >= d; i--)
    {
        /* p[i-d] = p[i-d] + c * p[i] */
        val0 = p[i-d];
        val1 = p[i];
        if (val0 >= n2)
            val0 -= n2;              /* [0, 2n) */
        umul_ppmm(p_hi, p_lo, c_precomp, val1);
        val1 = c * val1 - p_hi * n;  /* [0, 2n) */
        p[i-d] = val0 + val1;        /* [0, 4n) */
    }
}

/* division by x**d - 1 (not lazy: [0, n) -> [0, n)) */
FLINT_FORCE_INLINE
void _nmod_poly_divrem_circulant1(nn_ptr p, slong len, ulong d, ulong n)
{
    /* assumes len >= d */
    slong i;
    ulong j, r;

    r = len % d;
    i = len - r - d;  /* multiple of d, >= 0 by assumption */

    for (j = 0; j < r; j++)
        p[i+j] = n_addmod(p[i+j], p[d+i+j], n);

    i -= d;
    while (i >= 0)
    {
        for (j = 0; j < d; j++)
            p[i+j] = n_addmod(p[i+j], p[d+i+j], n);
        i -= d;
    }
}

/* assumes len > 0 and d > 0, multiples of 4 */
FLINT_FORCE_INLINE
void _nmod_poly_divrem_circulant1_v1(nn_ptr p, slong len, ulong d, ulong n)
{
    ulong i, j;

    for (j = d; j+d-1 < (ulong)len; j+=d)
    {
        for (i = 0; i < d; i+=4)
        {
            p[i+0] = n_addmod(p[i+0], p[j+i+0], n);
            p[i+1] = n_addmod(p[i+1], p[j+i+1], n);
            p[i+2] = n_addmod(p[i+2], p[j+i+2], n);
            p[i+3] = n_addmod(p[i+3], p[j+i+3], n);
        }
    }
    for (i = 0; i+j < (ulong)len; i+=4)
    {
        p[i+0] = n_addmod(p[i+0], p[j+i+0], n);
        p[i+1] = n_addmod(p[i+1], p[j+i+1], n);
        p[i+2] = n_addmod(p[i+2], p[j+i+2], n);
        p[i+3] = n_addmod(p[i+3], p[j+i+3], n);
    }
}

/*---------------------------------------------------*/
/* "c-circulant" division with remainder, transposed */
/*   -> expand sequence mod x**d - c                 */
/*---------------------------------------------------*/

/* transposed version of above function: 
 *    Input:
 *        :p: vector of length >= len
 *        :len: target length of expansion
 *        :d: positive integer
 *        :c: element of base field
 *    Effect:
 *    find coefficients d...len-1 of p by unrolling the recurrence with
 *    charpoly x**d - c
 *    -> explicitly: for all d*i+j with j < d such that d*i+j < len,
 *        p_{d*i+j} == p_j * c**i
 *   (in particular, first d entries are unchanged)
 **/
/* FIXME be more clear on what laziness is needed
 * in can be whatever
 * out is [0..2n) for new values; old values unchanged (needs reduction?)  */
void _nmod_poly_divrem_circulant_lazy_4_2_t(nn_ptr p, ulong len, ulong d, ulong c, ulong c_precomp, ulong n)
{
    ulong i, val, p_hi, p_lo;

    for (i = 0; i+d < len; i++)
    {
        /* p[i+d] = c * p[i] */
        val = p[i];
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        p[i+d] = c * val - p_hi * n;  /* [0, 2n) */
    }
}

void _nmod_poly_divrem_circulant1_t(nn_ptr p, ulong len, ulong d)
{
    ulong i;
    for (i = 0; i+d < len; i++)
        p[i+d] = p[i];
}


/*---------------------------------------------------*/
/* division with remainder mod product */
/*   TODO                 */
/*---------------------------------------------------*/

/*---------------------------------------------------*/
/* division with remainder mod product, transposed */
/*---------------------------------------------------*/

/** Division mod product of x**d - w, transposed
 *     Input:
 *         :p: array
 *         :len: positive integer (length of p)
 *         :d: positive integer
 *         :depth: nonnegative integer, current depth in root tree
 *         :node: nonnegative integer, current node in root tree
 *         :F: n_fft_args_t
 *     Requirements: d <= 2**depth, number of roots in F->tab_w at least 2**depth * node + d
 *     Action:
 *         store in p[d:len] the coefficients obtained by unrolling the
 *         length-d recurrence provided by
 *            F->prod(x - tab_w[2 * (node * 2**depth + k)] for k in range(d))
 *         on the initial d coefficients p[:d]
 *         (note that the input coefficients p[d:len] are ignored and overwritten)
 *     Algorithm:
 *     one could derive explanations similar to the comments at the beginning
 *     of TODO (name) reduce_mod_prod, or simply see it as its direct transpose
 */

void _nmod_poly_divrem_prod_roots_unity_t_lazy_x_x(nn_ptr p, ulong len, ulong d,
                                                   ulong depth, ulong node, n_fft_args_t F)
{
    /* base case: if len <= d, all requested terms are known */
    if (len <= d)
        return;

    /* currently, d <= 2**depth --> ensure 2**(depth-1) < d <= 2**depth */
    ulong depth_d = (d == 1) ? 0 : n_clog2_ge2(d);
    node = node << (depth - depth_d);
    depth = depth_d;
    ulong e = UWORD(1) << depth;

    /* if d is a power of 2, i.e. d == e (which covers the case d == 1 <=> depth == 0), */
    /* we are just unrolling with charpoly x**d - tab_w[2*node] */
    if (d == e)
    {
        _nmod_poly_divrem_circulant_lazy_4_2_t(p, len, d, F->tab_w[2*node], F->tab_w[2*node+1], F->mod);
        return;
    }

    ulong llen = FLINT_MIN(e, len);
    e = e/2;
    /* from here on, 1 <= e == 2**(depth-1) < d < llen <= 2**depth */

    const ulong w = F->tab_w[4*node];
    const ulong wpre = F->tab_w[4*node+1];
    ulong val0, val1, p_hi, p_lo;
    for (ulong i = 0; i < d - e; i++)
    {
        /* p[e + i] -= w * p[i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        umul_ppmm(p_hi, p_lo, wpre, val0);
        val0 = w * val0 - p_hi * F->mod;  /* [0, 2n) */
        p[e+i] = val1 + F->mod2 - val0;   /* [0, 4n) */
    }

    _nmod_poly_divrem_prod_roots_unity_t_lazy_x_x(p + e, llen-e, d-e, depth-1, 2*node+1, F);

    for (ulong i = 0; i < llen - e; i++)
    {
        /* p[e + i] += w * p[i] */
        val0 = p[i];
        val1 = p[e+i];
        if (val1 >= F->mod2)
            val1 -= F->mod2;              /* [0, 2n) */
        umul_ppmm(p_hi, p_lo, wpre, val0);
        val0 = w * val0 - p_hi * F->mod;  /* [0, 2n) */
        p[e+i] = val0 + val1;             /* [0, 4n) */
    }

    if (len > 2*e)
        _nmod_poly_divrem_circulant_lazy_4_2_t(p, len, 2*e, F->tab_w[2*node], F->tab_w[2*node+1], F->mod);
}


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

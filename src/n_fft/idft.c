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
#include "basic.c"

/*---------*/
/* helpers */
/*---------*/

// FIXME repeated from dft.c, see about making common basic macros / defs file
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

/*--------------*/
/* 2-point IDFT */
/*--------------*/

/** Gentleman-Sande butterfly, general
 * * in [0..2n) / out [0..2n) / max < 4n
 * * In-place transform
 *                            [1  w]
 *           [a  b] <- [a  b] [1 -w]
 * * n2 is 2*n, iw_pr is the precomputed data for multiplication by iw mod n
 *   p_hi, p_lo, tmp are temporaries
 * * can be seen as interpolation at points w = 1 / iw and -w, up to a scaling
 * by 1/2, since the inverse of [1  iw] is  1/2 * [1  1]
 *                              [1 -iw]           [w -w]
 */
#define IDFT2_LAZY22(a, b, n, n2, w, w_pr, p_hi, p_lo, tmp)   \
do {                                                          \
    tmp = (a) + (n2) - (b);     /* [0..4n) */                 \
    (a) = (a) + (b);            /* [0..4n) */                 \
    if ((a) >= (n2))                                          \
        (a) -= (n2);            /* [0..2n) */                 \
    N_MULMOD_PRECOMP_LAZY((b), w, tmp, w_pr, n, p_hi, p_lo);  \
                                /* --> (b) in [0..2n) */      \
} while(0)

/*--------------*/
/* 4-point IDFT */
/*--------------*/

/** 4-point IDFT, general
 * * in [0..2n) / out [0..2n) / max < 4n
 * * In-place transform
 *                              [1  w2    w1    w3]
 *                              [1 -w2  I*w1 -I*w3]
 * [a  b  c  d] <- [a  b  c  d] [1  w2   -w1   -w3]
 *                              [1 -w2 -I*w1  I*w3]
 *
 *                              [1     1   ] [1  w2        ]
 *                              [   1     I] [1 -w2        ]
 *              == [a  b  c  d] [1    -1   ] [       w1  w3]
 *                              [   1    -I] [       w1 -w3]
 */
#define IDFT4_LAZY22(a,b,c,d,                                      \
                     I,I_pr,w1,w1_pr,w2,w2_pr,w3,w3_pr,            \
                     n,n2,n4,p_hi,p_lo)                            \
do {                                                               \
        const ulong u0 = (a);                                      \
        const ulong u1 = (b);                                      \
        const ulong u2 = (c);                                      \
        const ulong u3 = (d);                                      \
                                                                   \
        ulong u4 = u0 + u2;            /* [0..4n) */               \
        ulong u5 = u0 + n2 - u2;       /* [0..4n) */               \
        ulong u6 = u1 + u3;            /* [0..4n) */               \
        ulong u7 = u1 + n2 - u3;       /* [0..4n) */               \
                                                                   \
        N_MULMOD_PRECOMP_LAZY(u7, I, u7, I_pr, n, p_hi, p_lo);     \
                                                                   \
        p_lo = u4 + u6;                /* [0..8n) */               \
        if (p_lo >= n4)                                            \
            p_lo -= n4;                                            \
        if (p_lo >= n2)                                            \
            p_lo -= n2;                                            \
        (a) = p_lo;                    /* [0..2n) */               \
                                                                   \
        u4 = u4 + n4 - u6;                                         \
        N_MULMOD_PRECOMP_LAZY((b), w2, u4, w2_pr, n, p_hi, p_lo);  \
        u6 = u5 + u7;                                              \
        N_MULMOD_PRECOMP_LAZY((c), w1, u6, w1_pr, n, p_hi, p_lo);  \
        u5 = u5 + n2 - u7;                                         \
        N_MULMOD_PRECOMP_LAZY((d), w3, u5, w3_pr, n, p_hi, p_lo);  \
} while(0)

void idft_lazy22(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 1)
    {
        ulong p_hi, p_lo, tmp;
        IDFT2_LAZY22(p[0], p[1], F->mod, F->mod2, F->tab_w[2*node], F->tab_w[2*node+1], p_hi, p_lo, tmp);
    }
    else
    {
        const ulong len = UWORD(1) << depth;
        idft_lazy22(p, depth-1, 2*node, F);
        idft_lazy22(p+len/2, depth-1, 2*node+1, F);

        const ulong w = F->tab_w[4*node];
        const ulong w_pr = F->tab_w[4*node+1];
        ulong p_hi, p_lo, tmp;

        for (ulong k = 0; k < len/2; k++)
        {
            IDFT2_LAZY22(p[k], p[len/2 + k], F->mod, F->mod2, w, w_pr, p_hi, p_lo, tmp);
        }
    }
}

void idft_node0_lazy12(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 0)
        return;

    if (depth == 1)
    {
        ulong tmp;
        DFT2_NODE0_LAZY12(p[0], p[1], F->mod, tmp);
    }
    else
    {
        const ulong len = UWORD(1) << depth;
        idft_node0_lazy12(p, depth-1, F);
        idft_lazy22(p+len/2, depth-1, 1, F);

        const ulong I = F->tab_w[0];
        const ulong I_pr = F->tab_w[1];
        ulong p_hi, p_lo, tmp;

        for (ulong k = 0; k < len/2; k++)
        {
            IDFT2_LAZY22(p[k], p[len/2 + k], F->mod, F->mod2, I, I_pr, p_hi, p_lo, tmp);
        }
    }
}

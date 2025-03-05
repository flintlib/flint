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

/*--------------*/
/* 2-point IDFT */
/*--------------*/

/** Gentleman-Sande butterfly, general
 * * in [0..2n) / out [0..2n) / max < 4n
 * * In-place transform
 *                            [1  iw]
 *           [a  b] <- [a  b] [1 -iw]
 * * n2 is 2*n, iw_pr is the precomputed data for multiplication by iw mod n
 *   p_hi, p_lo, tmp are temporaries
 * * can be seen as interpolation at points w = 1 / iw and -w, up to a scaling
 * by 1/2, since the inverse of [1  iw] is  1/2 * [1  1]
 *                              [1 -iw]           [w -w]
 */
// TODO make order of arguments consistent
#define IDFT2_LAZY22(a, b, n, n2, w, w_pr, p_hi, p_lo, tmp)   \
do {                                                          \
    tmp = (a) + (n2) - (b);     /* [0..4n) */                 \
    (a) = (a) + (b);            /* [0..4n) */                 \
    if ((a) >= (n2))                                          \
        (a) -= (n2);            /* [0..2n) */                 \
    N_MULMOD_PRECOMP_LAZY((b), w, tmp, w_pr, n, p_hi, p_lo);  \
                                /* --> (b) in [0..2n) */      \
} while(0)

// move in macros?
// in [0..4n) x [0..2n) -> out [0..4n) x [0..4n)
// TODO rename
#define BUTTERFLY_LAZY4244(a, b, n2, tmp)                     \
do {                                                          \
    tmp = (a);                                                \
    if (tmp >= (n2))                                          \
        tmp -= (n2);            /* [0..2n) */                 \
    (a) = tmp + (b);            /* [0..4n) */                 \
    (b) = tmp + (n2) - (b);     /* [0..4n) */                 \
} while(0)

/*--------------*/
/* 4-point IDFT */
/*--------------*/

/** 4-point IDFT, general
 * * in [0..4n) / out [0..4n) / max < 4n
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
 */
#define IDFT4_LAZY22(a, b, c, d,                                   \
                     w1, w1_pr, w2, w2_pr, w3, w3_pr,              \
                     n, n2, p_hi, p_lo)                            \
do {                                                               \
    const ulong v0 = (a);                                          \
    const ulong v1 = (b);                                          \
    const ulong v2 = (c);                                          \
    const ulong v3 = (d);                                          \
    ulong v4 = v0 + v1;                       /* < 4*n */          \
    if (v4 >= (n2))                                                \
        v4 -= (n2);                           /* < 2*n */          \
    ulong v5;                                                      \
    N_MULMOD_PRECOMP_LAZY(v5, (w2), v0 + (n2) - v1, (w2_pr), (n),  \
                            p_hi, p_lo);      /* < 2*n */          \
    ulong v6 = v2 + v3;                       /* < 4*n */          \
    if (v6 >= (n2))                                                \
        v6 -= (n2);                           /* < 2*n */          \
    ulong v7;                                                      \
    N_MULMOD_PRECOMP_LAZY(v7, (w3), v2 + (n2) - v3, (w3_pr), (n),  \
                            p_hi, p_lo);      /* < 2*n */          \
                                                                   \
    (a) = v4 + v6;                                                 \
    if ((a) >= (n2))                                               \
        (a) -= (n2);                           /* < 2*n */         \
    (b) = v5 + v7;                                                 \
    if ((b) >= (n2))                                               \
        (b) -= (n2);                           /* < 2*n */         \
    N_MULMOD_PRECOMP_LAZY((c), (w1), v4 + (n2) - v6, (w1_pr), (n), \
                            p_hi, p_lo);       /* < 2*n */         \
    N_MULMOD_PRECOMP_LAZY((d), (w1), v5 + (n2) - v7, (w1_pr), (n), \
                            p_hi, p_lo);       /* < 2*n */         \
} while(0)

#define IDFT4_LAZY12(a, b, c, d,                                   \
                     w1, w1_pr, w2, w2_pr, w3, w3_pr,              \
                     n, n2, p_hi, p_lo)                            \
do {                                                               \
    const ulong v0 = (a);                                          \
    const ulong v1 = (b);                                          \
    const ulong v2 = (c);                                          \
    const ulong v3 = (d);                                          \
    ulong v4 = v0 + v1;                       /* < 2*n */          \
    ulong v5;                                                      \
    N_MULMOD_PRECOMP_LAZY(v5, (w2), v0 + (n) - v1, (w2_pr), (n),   \
                            p_hi, p_lo);      /* < 2*n */          \
    ulong v6 = v2 + v3;                       /* < 2*n */          \
    ulong v7;                                                      \
    N_MULMOD_PRECOMP_LAZY(v7, (w3), v2 + (n) - v3, (w3_pr), (n),   \
                            p_hi, p_lo);      /* < 2*n */          \
                                                                   \
    (a) = v4 + v6;                            /* < 4*n */          \
    if ((a) >= (n2))                                               \
        (a) -= (n2);                           /* < 2*n */         \
    (b) = v5 + v7;                            /* < 4*n */          \
    if ((b) >= (n2))                                               \
        (b) -= (n2);                           /* < 2*n */         \
    N_MULMOD_PRECOMP_LAZY((c), (w1), v4 + (n2) - v6, (w1_pr), (n), \
                            p_hi, p_lo);       /* < 2*n */         \
    N_MULMOD_PRECOMP_LAZY((d), (w1), v5 + (n2) - v7, (w1_pr), (n), \
                            p_hi, p_lo);       /* < 2*n */         \
} while(0)

/*--------------*/
/* 8-point IDFT */
/*--------------*/

// TODO doc
#define IDFT8_NODE0_LAZY14(p0, p1, p2, p3, p4, p5, p6, p7,  \
                           mod, mod2, tab_w)                \
do {                                                        \
    ulong p_hi, p_lo, tmp;                                  \
                                                            \
    IDFT4_NODE0_LAZY14(p0, p1, p2, p3,                      \
                       tab_w[2], tab_w[3],                  \
                       mod, mod2, p_hi, p_lo);              \
    IDFT4_LAZY12(p4, p5, p6, p7,                            \
                tab_w[2], tab_w[3],                         \
                tab_w[4], tab_w[5],                         \
                tab_w[6], tab_w[7],                         \
                mod, mod2, p_hi, p_lo);                     \
                                                            \
    BUTTERFLY_LAZY4244(p0, p4, mod2, tmp);                  \
    BUTTERFLY_LAZY4244(p1, p5, mod2, tmp);                  \
    BUTTERFLY_LAZY4244(p2, p6, mod2, tmp);                  \
    BUTTERFLY_LAZY4244(p3, p7, mod2, tmp);                  \
} while(0)


/** 8-point IDFT
 * TODO clean, check laziness
 * * in [0..?n) / out [0..?n) / max < ?n
 */
#define DFT8_LAZY12(p0, p1, p2, p3, p4, p5, p6, p7,              \
                    node, mod, mod2, tab_w)                      \
do {                                                             \
    ulong p_hi, p_lo, tmp;                                       \
                                                                 \
    const ulong w = tab_w[2*(node)];                             \
    const ulong w_pr = tab_w[2*(node)+1];                        \
                                                                 \
    IDFT4_LAZY12(p0, p1, p2, p3,                                 \
                 tab_w[4*(node)], tab_w[4*(node)+1],             \
                 tab_w[8*(node)], tab_w[8*(node)+1],             \
                 tab_w[8*(node)+2], tab_w[8*(node)+3],           \
                 mod, mod2, p_hi, p_lo);                         \
                                                                 \
    IDFT4_LAZY12(p4, p5, p6, p7,                                 \
                 tab_w[4*(node)+2], tab_w[4*(node)+3],           \
                 tab_w[8*(node)+4], tab_w[8*(node)+5],           \
                 tab_w[8*(node)+6], tab_w[8*(node)+7],           \
                 mod, mod2, p_hi, p_lo);                         \
                                                                 \
    IDFT2_LAZY22(p0, p4, mod, mod2, w, w_pr, p_hi, p_lo, tmp);  \
    IDFT2_LAZY22(p1, p5, mod, mod2, w, w_pr, p_hi, p_lo, tmp);  \
    IDFT2_LAZY22(p2, p6, mod, mod2, w, w_pr, p_hi, p_lo, tmp);  \
    IDFT2_LAZY22(p3, p7, mod, mod2, w, w_pr, p_hi, p_lo, tmp);  \
} while(0)














/*--------------*/
/* general IDFT */
/*--------------*/


// TODO doc
// TODO make sure this is tested (code coverage: including for small depths)
void idft_lazy12(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 1)
    {
        ulong p_hi, p_lo, tmp;
        IDFT2_LAZY22(p[0], p[1], F->mod, F->mod2, F->tab_w[2*node], F->tab_w[2*node+1], p_hi, p_lo, tmp);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        IDFT4_LAZY12(p[0], p[1], p[2], p[3],
                     F->tab_w[2*node], F->tab_w[2*node+1],
                     F->tab_w[4*node], F->tab_w[4*node+1],
                     F->tab_w[4*node+2], F->tab_w[4*node+3],
                     F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 3)
    {
        DFT8_LAZY12(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                    node, F->mod, F->mod2, F->tab_w);
    }
    else
    {
        const ulong len = UWORD(1) << depth;

        // 4 recursive calls with depth-2
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/4;
        const nn_ptr p2 = p + 2*len/4;
        const nn_ptr p3 = p + 3*len/4;
        idft_lazy12(p0, depth-2, 4*node, F);
        idft_lazy12(p1, depth-2, 4*node+1, F);
        idft_lazy12(p2, depth-2, 4*node+2, F);
        idft_lazy12(p3, depth-2, 4*node+3, F);

        const ulong w2 = F->tab_w[2*node];
        const ulong w2_pr = F->tab_w[2*node+1];
        const ulong w = F->tab_w[4*node];
        const ulong w_pr = F->tab_w[4*node+1];
        const ulong Iw = F->tab_w[4*node+2];
        const ulong Iw_pr = F->tab_w[4*node+3];
        ulong p_hi, p_lo;

        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_LAZY22(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY22(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY22(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY22(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
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
        BUTTERFLY_LAZY12(p[0], p[1], F->mod, tmp);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        IDFT4_NODE0_LAZY14(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3],
                           F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 3)
    {
        // TODO to be improved
        IDFT8_NODE0_LAZY14(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                           F->mod, F->mod2, F->tab_w);
    }
    else
    {
        const ulong len = UWORD(1) << depth;

        // 4 recursive calls with depth-2
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/4;
        const nn_ptr p2 = p + 2*len/4;
        const nn_ptr p3 = p + 3*len/4;
        idft_node0_lazy12(p0, depth-2, F);
        idft_lazy12(p1, depth-2, 1, F);
        idft_lazy12(p2, depth-2, 2, F);
        idft_lazy12(p3, depth-2, 3, F);

        // 4-point butterflies
        // input p0 in [0,4n), p1,p2,p3 in [0,2n)
        // output p0,p1,p2,p3 in [0,4n)
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_NODE0_LAZY4222(p0[k+0], p1[k+0], p2[k+0], p3[k+0],
                                 F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE0_LAZY4222(p0[k+1], p1[k+1], p2[k+1], p3[k+1],
                                 F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE0_LAZY4222(p0[k+2], p1[k+2], p2[k+2], p3[k+2],
                                 F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE0_LAZY4222(p0[k+3], p1[k+3], p2[k+3], p3[k+3],
                                 F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
        }
    }
}

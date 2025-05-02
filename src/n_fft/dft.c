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

/** Structure:
 * - the core function is dft_node_lazy_4_4, which goes down the subproduct
 *   tree from an arbitrary node in this tree; it takes input values in [0..4n)
 *   and return values in [0..4n)
 * TODO add/improve explanations below
 * - this core function costs more than a DFT at node 0 (at least for small /
 *   small-ish lengths), so a specific function for that are given, targeting
 *   input values in [0..n) and return values in [0..4n)
 *   (and in [0..2n) respectively (the former
 *   calls the latter)
 * - less lazy variants (in the end, this is often used
 *   non-lazy, [0..n) -> [0..n))
 */

/** 2**depth-point DFT, general node
 * * In-place transform p of length len == 2**depth, seen as a polynomial of
 * degree < len, into the concatenation of all polynomial evaluations
 *          [p(w_k), p(-w_k)] for k in range(len),
 * where w_k = F->tab_w[2**depth * node + 2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the len roots of the
 * polynomial x**len - F->tab_w[node]
 * * Requirement (not checked):
 *        3 <= depth
 *        (node+1) * 2**depth <= 2**F.depth (length of F->tab_w)
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
void dft_node_lazy_4_4(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 3)
    {
        DFT8_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        DFT16_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                            p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                            node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
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
            DFT4_NODE_LAZY_4_4(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_NODE_LAZY_4_4(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_NODE_LAZY_4_4(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
            DFT4_NODE_LAZY_4_4(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2, p_hi, p_lo, tmp);
        }

        // 4 recursive calls with depth-2
        dft_node_lazy_4_4(p0, depth-2, 4*node, F);
        dft_node_lazy_4_4(p1, depth-2, 4*node+1, F);
        dft_node_lazy_4_4(p2, depth-2, 4*node+2, F);
        dft_node_lazy_4_4(p3, depth-2, 4*node+3, F);
    }
}

/** 2**depth-point DFT
 * * In-place transform p of length len == 2**depth, seen as a polynomial of
 * degree < len, into the concatenation of all polynomial evaluations
 *          [p(w_k), p(-w_k)] for k in range(len),
 * where w_k = F->tab_w[2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the roots of the polynomial
 * x**len - 1, precisely they are all powers of the chosen len-th primitive
 * root of unity with exponents listed in bit reversed order
 * * Requirement (not checked): 3 <= depth <= F.depth
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 */
void dft_lazy_2_4(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 3)
    {
        DFT8_LAZY_2_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        DFT16_LAZY_2_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                       p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                       F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_LAZY_2_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
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
            DFT4_LAZY_2_4(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            if (p0[k] >= F->mod2)
                p0[k] -= F->mod2;
        }

        // 4 recursive calls with depth-2
        dft_lazy_2_4(p0, depth-2, F);
        dft_node_lazy_4_4(p1, depth-2, 1, F);
        dft_node_lazy_4_4(p2, depth-2, 2, F);
        dft_node_lazy_4_4(p3, depth-2, 3, F);
    }
}

void dft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 4)
    {
        DFT16_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                       p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                       F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        DFT32_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
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
            DFT4_LAZY_1_4(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            if (p0[k] >= F->mod2)
                p0[k] -= F->mod2;
        }

        // 4 recursive calls with depth-2
        dft_lazy_2_4(p0, depth-2, F);
        dft_node_lazy_4_4(p1, depth-2, 1, F);
        dft_node_lazy_4_4(p2, depth-2, 2, F);
        dft_node_lazy_4_4(p3, depth-2, 3, F);
    }
    else if (depth == 3)
    {
        DFT8_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        DFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
    }
    else if (depth == 1)
    {
        ulong tmp;
        DFT2_LAZY_1_2(p[0], p[1], F->mod, tmp);
    }
}

/*---------------*/
/* some comments */
/*---------------*/

/** Lazier variants for DFT with general node:
 * - lazy_1_4 variants would be basically identical to the lazy_2_4 variants (see the macros)
 * - writing lazy_2_4 variants of the DFTxx_NODE_LAZY_4_4 macros and then of
 * dft_node_lazy_4_4 brings almost no speedup (very marginal gain up to length
 * 32 or 64, nothing observable beyond this)
 */

/** Base cases:
 * - having macros for "small" lengths (up to 16 or 32 at least) improves performance
 * - removing the base cases depth==3 in internal functions where this case is
 *   not really used (eg dft_node_lazy_4_4) does not make a difference
 */

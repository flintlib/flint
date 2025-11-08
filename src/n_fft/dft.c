/*
    Copyright (C) 2024, 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft.h"
#include "n_fft/impl.h"
#include "n_fft/impl_macros_dft.h"
#include "n_fft/impl_macros_tft.h"

/** Structure for DFT.
 * 
 * The DFT functions take as input an array p of length ``len == 2**depth``,
 * and in-place modify it into the corresponding ``len`` DFT evaluations (see
 * functions documentations below for more details).
 *
 * - The main interface is n_fft_dft, it solves the problem at node 0
 *   (evaluating at all roots of unity of order 2**depth), as documented
 *   in n_fft.h.
 * - The core function is `dft_node_lazy_4_4`, which goes down the subproduct
 *   tree from an arbitrary node in this tree; it takes input values in [0..4n)
 *   and return values in [0..4n), following the idea of lazy butterflies
 *   highlighted by David Harvey [Faster arithmetic for number-theoretic
 *   transforms, Journal of Symbolic Computation, Volume 60, 2014, pp 113-119].
 * - This core function costs more than a DFT at node 0, at least for small or
 *   smallish lengths. So a specific function for node 0 is given
 *   (`dft_lazy_1_4`), targeting input values in [0..n) and return values in
 *   [0..4n) (it iself uses a similar function `dft_lazy_2_4`). The main
 *   function `n_fft_dft` just calls `dft_lazy_1_4` and then reduces the output
 *   to [0..n).
 */

/** Structure for TFT (very similar).
 *
 * Truncated variants of the DFT provide finer control and performance
 * depending on the input polynomial length (ilen) and output number of
 * evaluations (olen).
 *
 * The TFT functions below take as input an array of length max(ilen, len),
 *    p = [p_0,...,p_{ilen-1},xx,...,xx]
 * where len is the smallest power of 2 which is >= olen. Here each `xx` can be
 * any value, these values fill the array up to the required length and might
 * be modified during the algorithm.
 *
 * In output, the first olen entries of p contain the "first" DFT evaluations
 * (see functions documentations below for more details).
 *
 * - The main interface is n_fft_tft, it solves the problem at node 0
 *   (evaluating at roots of unity of order 2**depth), as documented
 *   in n_fft.h.
 * - The core function is `tft_node_lazy_4_4`, which goes down the subproduct
 *   tree from an arbitrary node in this tree; it takes input values in [0..4n)
 *   and return values in [0..4n), following the idea of lazy butterflies
 *   highlighted by David Harvey [Faster arithmetic for number-theoretic
 *   transforms, Journal of Symbolic Computation, Volume 60, 2014, pp 113-119].
 * - This core function costs more than a TFT at node 0, at least for small or
 *   smallish lengths. So a specific function for node 0 is given
 *   (`tft_lazy_1_4`), targeting input values in [0..n) and return values in
 *   [0..4n). The main function `n_fft_tft` just calls `tft_lazy_1_4` and then
 *   reduces the output to [0..n).
 */

/** Example for nodes/depth:
 *   if F->depth is 3, the tree of roots of unity in F->tab_w is
 *                    1                               d3n0                <-- depth 3 
 *               /        \                        /        \
 *             1            -1                 d2n0          d2n1         <-- depth 2
 *           /   \        /     \     =        /   \        /     \
 *         1     -1      I      -I         d1n0   d1n1   d1n2    d1n3     <-- depth 1
 *        / \    / \    / \    /  \         / \    / \    / \    /  \
 *       1  -1  I  -I  J  -J  IJ -IJ       1  -1  I  -I  J  -J  IJ -IJ    <-- depth 0
 *  stored as, ommitting precomputations:
 *    F->tab_w == [1, 1_pr, I, I_pr, J, J_pr, IJ, IJ_pr]
 *  (the elements -1, -I, -J, -IJ are not stored)
 *
 *
 *  -> calling a function with depth==3 and node==0 is performing
 *  evaluation at all these 8 points (8th roots of 1)
 *  -> calling a function with depth==2 and node==0 is performing
 *  evaluation at all points at the leaves of the left child d2n0
 *  of the root of the tree d3n0 (4th roots of 1)
 *  -> calling a function with depth==2 and node==1 is performing
 *  evaluation at all points at the leaves of the right child d2n1
 *  of d3n0 (4th roots of -1)
 *  -> calling a function with depth==1 and node==1 is performing
 *  evaluation at all points at the leaves of the subtree rooted
 *  at d1n1 (square roots of -1)
 *  -> calling a function with depth==1 and node==2 is performing
 *  evaluation at all points at the leaves of the subtree rooted
 *  at d1n2 (square roots of I)
 */

/*----------------------------*/
/*  DFT: auxiliary functions  */
/*----------------------------*/

/** 2**depth-point DFT, general node
 * * In-place transform p of length len == 2**depth, seen as a polynomial of
 * degree < len, into the concatenation of all polynomial evaluations
 *          [p(w_k), p(-w_k)] for k in range(len/2),
 * where w_k = F->tab_w[2**depth * node + 2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the len roots of the
 * polynomial x**len - F->tab_w[2*node] (for example, if depth=
 * * Requirements (not checked):
 *        3 <= depth
 *        (node+1) * 2**depth < 2**F->depth (length of F->tab_w)
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
void dft_node_lazy_4_4(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 2)  /* TODO added base case, check it does not impact performance */
    {
        DFT4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3],
                           F->tab_w[2*node+0], F->tab_w[2*node+1],    
                           F->tab_w[4*node+0], F->tab_w[4*node+1],   
                           F->tab_w[4*node+2], F->tab_w[4*node+3],
                           F->mod, F->mod2);
    }
    else if (depth == 3)
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

        for (ulong k = 0; k < len/4; k+=4)
        {
            DFT4_NODE_LAZY_4_4(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2);
            DFT4_NODE_LAZY_4_4(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2);
            DFT4_NODE_LAZY_4_4(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2);
            DFT4_NODE_LAZY_4_4(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2pre, w, wpre, Iw, Iwpre, F->mod, F->mod2);
        }

        // 4 recursive calls with depth-2
        dft_node_lazy_4_4(p0, depth-2, 4*node, F);
        dft_node_lazy_4_4(p1, depth-2, 4*node+1, F);
        dft_node_lazy_4_4(p2, depth-2, 4*node+2, F);
        dft_node_lazy_4_4(p3, depth-2, 4*node+3, F);
    }
}

/** 2**depth-point DFT
 * Same specification as n_fft_dft, except for:
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 *             requirement (not checked): depth <= F->depth
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 *             requirement (not checked): 3 <= depth <= F->depth
 */
void dft_lazy_2_4(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 2)  /* TODO added base case, check it does not impact performance */
    {
        DFT4_LAZY_2_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
    }
    else if (depth == 3)
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
        for (ulong k = 0; k < len/4; k++)
        {
            DFT4_LAZY_2_4(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
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
        for (ulong k = 0; k < len/4; k++)
        {
            DFT4_LAZY_1_4(p0[k], p1[k], p2[k], p3[k], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
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
        DFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
    }
    else if (depth == 1)
    {
        DFT2_LAZY_1_2(p[0], p[1], F->mod);
    }
}

/*----------------------------*/
/*  TFT: auxiliary functions  */
/*----------------------------*/

/* TODO fewer base cases in node_lazy_4_4, more base cases in lazy_1_4 */

/** truncated Fourier transform, ilen == len, general node
 * * same as dft_node_lazy_4_4, also with its requirements on depth and node,
 * but only computes the first olen evaluations
 * * Input requirements:
 *        ilen is len == 2**depth, power of 2 at least 4
 *        0 < olen <= len is a multiple of 4
 * * lazy_4_4: in [0..4n) / out [0..4n) / max < 4n
 */
void tft_node_lazy_4_4(nn_ptr p, ulong olen, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 2)
    {
        if (olen == 4)
        {
            DFT4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3],
                               F->tab_w[2*node], F->tab_w[2*node+1],
                               F->tab_w[4*node], F->tab_w[4*node+1],
                               F->tab_w[4*node+2], F->tab_w[4*node+3],
                               F->mod, F->mod2);
        }
        else
        {
            TFT4_2_NODE_LAZY_4_4(p[0], p[1], p[2], p[3],
                                 F->tab_w[2*node], F->tab_w[2*node+1],
                                 F->tab_w[4*node], F->tab_w[4*node+1],
                                 F->mod, F->mod2);
        }
    }
    else if (depth == 3)
    {
        if (olen == 8)
        {
            DFT8_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                               node, F->mod, F->mod2, F->tab_w);
        }
        else
        {
            TFT8_4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                 node, F->mod, F->mod2, F->tab_w);
        }
    }
    else if (depth == 4)
    {
        if (olen == 16)
        {
            DFT16_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 4)
        {
            TFT16_4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                  p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                  node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 8)
        {
            TFT16_8_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                  p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                  node, F->mod, F->mod2, F->tab_w);
        }
        else  /* olen == 12 */
        {
            TFT16_12_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   node, F->mod, F->mod2, F->tab_w);
        }
    }
    else if (depth == 5)
    {
        if (olen == 32)
        {
            DFT32_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 4)
        {
            TFT32_4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                  p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                  p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                  p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                  node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 8)
        {
            TFT32_8_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                  p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                  p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                  p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                  node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 12)
        {
            TFT32_12_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                   p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                   node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 16)
        {
            TFT32_16_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                   p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                   node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 20)
        {
            TFT32_20_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                   p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                   node, F->mod, F->mod2, F->tab_w);
        }
        else if (olen == 24)
        {
            TFT32_24_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                   p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                   node, F->mod, F->mod2, F->tab_w);
        }
        else  /* olen == 28 */
        {
            TFT32_28_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                                   p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
                                   node, F->mod, F->mod2, F->tab_w);
        }
    }
    else if (olen == (UWORD(1) << depth))
    {
        dft_node_lazy_4_4(p, depth, node, F);
    }
    else if (olen <= (UWORD(1) << (depth - 1)))
    {
        const ulong len = UWORD(1) << depth;
        ulong new_depth = n_clog2_ge2(olen);
        node = node << (depth - new_depth);
        depth = new_depth;

        /* reduce p1 mod x**(len_rec) - root */
        ulong len_rec = UWORD(1) << depth;
        _nmod_poly_divrem_circulant_lazy_4_4(p, len, len_rec, F->tab_w[node], F->tab_w[node+1], F->mod, F->mod2);
        tft_node_lazy_4_4(p, olen, depth, node, F);
    }
    else
    {
        /* current full length */
        /* -> input assumption: len/2 < olen < len */
        const ulong len = UWORD(1) << depth;

        /* two recursive calls in len/2 */
        // in: [0..4n), out: [0..4n)
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        const ulong w = F->tab_w[2*node];
        const ulong wpre = F->tab_w[2*node+1];
        for (ulong k = 0; k < len/2; k+=4)
        {
            DFT2_NODE_LAZY_4_4(p0[k+0], p1[k+0], w, wpre, F->mod, F->mod2);
            DFT2_NODE_LAZY_4_4(p0[k+1], p1[k+1], w, wpre, F->mod, F->mod2);
            DFT2_NODE_LAZY_4_4(p0[k+2], p1[k+2], w, wpre, F->mod, F->mod2);
            DFT2_NODE_LAZY_4_4(p0[k+3], p1[k+3], w, wpre, F->mod, F->mod2);
        }

        dft_node_lazy_4_4(p0, depth-1, 2*node, F);
        tft_node_lazy_4_4(p1, olen - len/2, depth-1, 2*node+1, F);
    }
}

/** truncated Fourier transform, general ilen, node 0
 * * similar to dft_lazy_1_4, but accepts p of arbitrary length and only
 * computes the first olen evaluations
 * * same requirements on depth and node as in dft_lazy_1_4, where here depth
 * is understood as ceiling(log_2(olen))
 * * Input requirements:
 *        ilen is a multiple of 4, at least 8
 *        olen is a positive multiple of 4
 *           (exception: ilen==4 is accepted when olen==4)
 *        p has space for max(ilen, len) coefficients where len = 1<<depth, and
 *        if ilen < len then coefficients ilen..len-1 may be modified by the
 *        algorithm
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 *
 * NOTE: some calls could be lazier, but the gain should be quite modest:
 *     [comment1]   dft_node_lazy_4_4 -> could benefit from lazy_2_4
 *     [comment2]  DFT2_NODE_LAZY_4_4 -> could benefit from LAZY_1_4
 */
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F)
{
    /* some recursive calls arrive at this base case */
    if (olen == 4 && ilen == 4)
    {
        DFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
        return;
    }

    const ulong odepth = n_clog2_ge2(olen);
    const ulong idepth = n_clog2_ge2(ilen);

    /* ilen > len: reduce mod x**len - 1 and call again */
    if (idepth > odepth)
    {
        const ulong len = UWORD(1) << odepth;
        _nmod_poly_divrem_circulant1(p, ilen, len, F->mod);
        tft_lazy_1_4(p, len, olen, F);
    }

    /* ilen ~ olen: recurse into dft + tft at half length */
    else if (idepth == odepth)
    {
        const ulong len = UWORD(1) << odepth;
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        ulong k;
        for (k = 0; k < ilen - len/2; k++)
            DFT2_LAZY_1_2(p0[k], p1[k], F->mod);
        for (; k < len/2; k++)
            p1[k] = p0[k];
        dft_lazy_2_4(p0, odepth - 1, F);
        tft_node_lazy_4_4(p1, olen - len/2, odepth - 1, 1, F);
    }

    /* if ilen ~ olen/2, do special butterfly and recurse */
    /* this is not necessary: it is covered by the next (convoluted) "if" */
    /* but since this case is frequent, let's write it in a simple way */
    else if (idepth == odepth - 1)
    {
        const ulong ilen2 = UWORD(1) << idepth;
        const ulong len_rec = ilen2 / 2;
        if (olen <= ilen2 + len_rec)
        {
            /* olen in [ilen2 + 4, ilen2 + len_rec] */
            /* -> a single TFT call */
            const nn_ptr p0 = p + 0*len_rec;
            const nn_ptr p1 = p + 1*len_rec;
            const nn_ptr p2 = p + 2*len_rec;
            ulong k, tmp0, tmp1;

            /* butterflies */
            for (k = 0; k < ilen - len_rec; k++)
            {
                tmp0 = p0[k];
                tmp1 = p1[k];
                DFT2_LAZY_1_2(p0[k], p1[k], F->mod);
                TFT2_1_NODE_LAZY_4_4(tmp0, tmp1, F->tab_w[2], F->tab_w[3], F->mod, F->mod2);  /* [comment2] */
                p2[k] = tmp0;
            }
            for (; k < len_rec; k++)
            {
                tmp0 = p0[k];
                p1[k] = tmp0;
                p2[k] = tmp0;
            }

            /* recursive calls */
            dft_lazy_2_4(p0, idepth - 1, F);
            dft_node_lazy_4_4(p1, idepth - 1, 1, F);
            tft_node_lazy_4_4(p2, olen - ilen2, idepth - 1, 2, F);
        }
        else
        {
            /* olen in [ilen2 + len_rec + 4, len] */
            /* -> one DFT call + one TFT call */
            const nn_ptr p0 = p + 0*len_rec;
            const nn_ptr p1 = p + 1*len_rec;
            const nn_ptr p2 = p + 2*len_rec;
            const nn_ptr p3 = p + 3*len_rec;
            ulong k, tmp0, tmp1;

            /* butterflies */
            for (k = 0; k < ilen - len_rec; k++)
            {
                tmp0 = p0[k];
                tmp1 = p1[k];
                DFT2_LAZY_1_2(p0[k], p1[k], F->mod);
                DFT2_NODE_LAZY_4_4(tmp0, tmp1, F->tab_w[2], F->tab_w[3], F->mod, F->mod2);  /* [comment2] */
                p2[k] = tmp0;
                p3[k] = tmp1;
            }
            for (; k < len_rec; k++)
            {
                /* butterfly with input p1[k] == 0: p0[k] unchanged, p1[k] <- p[k] */
                tmp0 = p0[k];
                p1[k] = tmp0;
                p2[k] = tmp0;
                p3[k] = tmp0;
            }
            dft_lazy_2_4(p0, idepth - 1, F);
            dft_node_lazy_4_4(p1, idepth - 1, 1, F);
            dft_node_lazy_4_4(p2, idepth - 1, 2, F);
            tft_node_lazy_4_4(p3, olen - ilen2 - len_rec, idepth - 1, 3, F);
        }
    }

    else  /* idepth < odepth - 1 */
    {
        const ulong ilen2 = UWORD(1) << idepth;
        ulong i, k;
        ulong val;
        nn_ptr pp = p + ilen2;
        const ulong len_rec = ilen2 / 2;

        /* first DFT round with node == 0 */
        /* save values and perform butterflies */
        for (k = 0; k < ilen - len_rec; k++)
        {
            pp[k] = p[k];
            pp[len_rec + k] = p[len_rec + k];
            DFT2_LAZY_1_2(p[k], p[len_rec + k], F->mod);
        }
        for (; k < len_rec; k++)
        {
            /* butterfly with input p[len_rec + k] == 0: p[k] unchanged, p[len_rec + k] = p[k] */
            val = p[k];
            p[len_rec + k] = val;
            pp[k] = val;
        }
        dft_lazy_2_4(p, idepth - 1, F);
        dft_node_lazy_4_4(p+len_rec, idepth - 1, 1, F);  /* [comment1] */

        /* update */
        p += ilen2;
        pp += ilen2;
        olen -= ilen2;

        /* next DFT rounds */
        for (i = 1; olen > ilen2; i++)
        {
            /* save values and perform butterflies */
            for (k = 0; k < ilen - len_rec; k++)
            {
                pp[k] = p[k];
                pp[len_rec + k] = p[len_rec + k];
                DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment2] */
                                   F->tab_w[2*i], F->tab_w[2*i+1],
                                   F->mod, F->mod2);
            }
            for (; k < len_rec; k++)
            {
                /* butterfly with input p[len_rec + k] == 0: p[k] unchanged, p[len_rec + k] = p[k] */
                val = p[k];
                p[len_rec + k] = val;
                pp[k] = val;
            }
            dft_node_lazy_4_4(p, idepth - 1, 2*i, F);
            dft_node_lazy_4_4(p+len_rec, idepth - 1, 2*i+1, F);

            p += ilen2;
            pp += ilen2;
            olen -= ilen2;
        }

        /* butterflies */
        if (olen > len_rec)
        {
            for (k = 0; k < ilen - len_rec; k++)
            {
                DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment2] */
                                F->tab_w[2*i], F->tab_w[2*i+1],
                                F->mod, F->mod2);
            }
            for (; k < len_rec; k++)
                p[len_rec + k] = p[k];
            dft_node_lazy_4_4(p, idepth - 1, 2*i, F);
            tft_node_lazy_4_4(p+len_rec, olen - len_rec, idepth - 1, 2*i+1, F);
        }
        else
        {
            for (k = 0; k < ilen - len_rec; k++)
            {
                TFT2_1_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment2] */
                                     F->tab_w[2*i], F->tab_w[2*i+1],
                                     F->mod, F->mod2);
            }
            tft_node_lazy_4_4(p, olen, idepth - 1, 2*i, F);
        }
    }
}

/*-------------------*/
/*  main interfaces  */
/*-------------------*/

/* depth > 0, p has length 1<<depth */
void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    n_fft_args_t Fargs;
    n_fft_set_args(Fargs, F->mod, F->tab_w);
    dft_lazy_1_4(p, depth, Fargs);
    for (ulong k = 0; k < (UWORD(1) << depth); k++)
    {
        if (p[k] >= Fargs->mod2)
            p[k] -= Fargs->mod2;
        if (p[k] >= Fargs->mod)
            p[k] -= Fargs->mod;
    }
}

/* depth > 0, p has length 1<<depth */
void n_fft_idft_t(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    n_fft_args_t Fargs;
    n_fft_set_args(Fargs, F->mod, F->tab_iw);
    dft_lazy_1_4(p, depth, Fargs);

    // see comments in idft concerning this loop
    const ulong inv2 = F->tab_inv2[2*depth-2];
    const ulong inv2_pr = F->tab_inv2[2*depth-1];
    for (ulong k = 0; k < (UWORD(1) << depth); k++)
        p[k] = n_mulmod_shoup(inv2, p[k], inv2_pr, F->mod);
}

/* olen > 0, ilen >= 0, both multiples of 4 */
/* p has length max(ilen, next_power_of_two(olen)) */
void n_fft_tft(nn_ptr p, ulong ilen, ulong olen, n_fft_ctx_t F)
{
    if (ilen == 0)
    {
        for (ulong k = 0; k < olen; k++)
            p[k] = 0;
        return;
    }

    n_fft_args_t Fargs;
    n_fft_set_args(Fargs, F->mod, F->tab_w);
    tft_lazy_1_4(p, ilen, olen, Fargs);
    for (ulong k = 0; k < olen; k++)
    {
        if (p[k] >= Fargs->mod2)
            p[k] -= Fargs->mod2;
        if (p[k] >= Fargs->mod)
            p[k] -= Fargs->mod;
    }
}

/*---------------*/
/* some comments */
/*---------------*/

/** In n_fft_idft_t, there is apparently no gain from using the lazy
 * mulmod_shoup variant whose output is in [0..2n) (so one may as well use the
 * non-lazy one which ensures output < n)              
 */

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

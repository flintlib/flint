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

/** Structure.
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

/** Example for nodes/depth:
 *   if F.depth is 3, the tree of roots of unity in F->tab_w is
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

/*-----------------------*/
/*  auxiliary functions  */
/*-----------------------*/

/** 2**depth-point DFT, general node
 * * In-place transform p of length len == 2**depth, seen as a polynomial of
 * degree < len, into the concatenation of all polynomial evaluations
 *          [p(w_k), p(-w_k)] for k in range(len),
 * where w_k = F->tab_w[2**depth * node + 2*k] for 0 <= k < 2**(depth-1)
 * * By construction these evaluation points are the len roots of the
 * polynomial x**len - F->tab_w[2*node] (for example, if depth=
 * * Requirements (not checked):
 *        3 <= depth
 *        (node+1) * 2**depth < 2**F.depth (length of F->tab_w)
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
 *             requirement (not checked): depth <= F.depth
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 *             requirement (not checked): 3 <= depth <= F.depth
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

/*-------------------*/
/*  main interfaces  */
/*-------------------*/

void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth > 0)
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
}

void n_fft_idft_t(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth > 0)
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

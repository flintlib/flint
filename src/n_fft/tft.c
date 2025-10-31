/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft.h"
#include "dft.c" // TODO put dft_node and others in n_fft.h
#include "n_fft_macros.h"
#include "nmod_poly.h"

/** Structure.
 * - The main interface is n_fft_tft, it solves the problem at node 0
 *   (evaluating at roots of unity of order 2**depth), as documented
 *   in n_fft.h.
 * - UPDATE The core function is `tft_node_lazy_4_4`, which goes down the subproduct
 *   tree from an arbitrary node in this tree; it takes input values in [0..4n)
 *   and return values in [0..4n), following the idea of lazy butterflies
 *   highlighted by David Harvey [Faster arithmetic for number-theoretic
 *   transforms, Journal of Symbolic Computation, Volume 60, 2014, pp 113-119].
 * - UPDATE This core function costs more than a DFT at node 0, at least for small or
 *   smallish lengths. So a specific function for node 0 is given
 *   (`dft_lazy_1_4`), targeting input values in [0..n) and return values in
 *   [0..4n) (it iself uses a similar function `dft_lazy_2_4`). The main
 *   function `n_fft_dft` just calls `dft_lazy_1_4` and then reduces the output
 *   to [0..n).
 */

/** For explanations about nodes/depth, see dft.c */

/** TODO explanations of ilen/olen */




/** truncated Fourier transform
Input:
    :f: coefficients [p_0,...,p_{ilen-1},xx,...,xx] of some polynomial of degree < ilen
            if ilen < NextPowerOfTwo(olen), xx,...,xx are any values to reach
            this power of 2 length
            (these values may be overwritten)
    :ilen: positive integer, strict degree bound for p
    :olen: positive integer, number of requested output values
    :depth: nonnegative integer, current depth in recursive tree
    :node: nonnegative integer, current root in recursive tree
Requirements:
    olen <= 2**depth
    (node+1) * 2**depth <= precomputation length
Output:
    in-place modify p so that it holds the first olen bit-reversed
        FFT evaluations p(w**0),w**(d/2),w**(d/4),w**(3d/4),...  (for node == 0);
    more precisely it holds the sums
        sum(p[k] * ws_br[node * 2**depth + i]**k for k in range(ilen)), for i in range(olen)
        --> TODO simply say that it changes first olen coeffs into the first olen ones as if running the full DFT, and refer to doc in dft.c ??
Algo:
    uses a direct reduction-tree approach:
    --> we are currently at node x**d - ws_br[node], where d = 2**depth
    --> subtrees have as roots the nodes x**(d/2) - ws_br[2*node] and x**(d/2) - ws_br[2*node+1]
**/
void n_fft_tft_draft(nn_ptr p, ulong ilen, ulong olen, ulong depth, ulong node, n_fft_args_t F)
{
    /* current full length */
    const ulong len = UWORD(1) << depth;

    /* if olen <= len/2, go down the tree */
    if (2*olen <= len)
        n_fft_tft_draft(p, ilen, olen, depth-1, 2*node, F);

    /* -> in what follows, len/2 < olen <= len */

    if (olen == 1)  /* also means len == 1: single point evaluation */
    {
        /* FIXME as such, this requires max(p) + 2n - 1 < 2**FLINT_BITS */
        /* FIXME efficiency: this evaluation should be done faster for point 1 or -1 */
        ulong val = _nmod_poly_evaluate_nmod_precomp_lazy(p, ilen, F->tab_w[2*node], F->tab_w[2*node+1], F->mod);
        return;
    }
    /* FIXME do olen == 2 similarly? or more base cases here for small olen? */
    /* or will none of these small olen be used anyway, due to some roundup? */

    if (ilen == 1)  /* this case is not necessary (the above one is) */
    {
        for (ulong k = 1; k < olen; k++)
            p[k] = p[0];
        return;
    }

    /* now len >= 2, use two recursive calls in length len/2 */
    const ulong len_rec = len >> 1;  // just use len/2 everywhere?
    const ulong olen_rec = olen - len_rec;

    /* if ilen <= len/2, directly do recursive calls */
    /* (this is a simplification of the general case below, and also */
    /* allows to do both calls with ilen again instead of len_rec) */
    /* TODO think about changing this into a loop beforehand,
     * which would ensures we have len/2 < ilen < len,
     * and then apply one butterfly which implies ilen == (o?)len?! */
    if (2*ilen <= len)
    {
        for (ulong k = 0; k < ilen; k++)
            p[len_rec+k] = p[k];
        n_fft_tft_draft(p        , ilen, len_rec , depth-1, 2*node  , F);
        n_fft_tft_draft(p+len_rec, ilen, olen_rec, depth-1, 2*node+1, F);
        return;
    }

    /* if ilen > len, reduce f mod x**len - root, and continue */
    if (ilen > len)
    {
        reduce_mod_xnma(p, ilen, len, F->tab_w[2*node]);
        ilen = len;
    }

    /* now len/2 < ilen <= len, do butterfly and then recursive calls */
    /* ww = ws_br[2*node] */
    /* for k in range(ilen - len_rec): */
    /*     tmp = ww * f[len_rec + k] */
    /*     f[len_rec + k] = f[k] - tmp */
    /*     f[k] = f[k] + tmp */
    /* for k in range(ilen - len_rec, len_rec): */
    /*     f[len_rec + k] = f[k] */

    // in: [0..4n), out: [0..4n)
    const nn_ptr p0 = p;
    const nn_ptr p1 = p + len_rec;
    const ulong w = F->tab_w[2*node];
    const ulong wpre = F->tab_w[2*node+1];
    for (ulong k = 0; k < len/2; k+=4)
    {
        DFT2_NODE_LAZY_4_4(p0[k+0], p1[k+0], w, wpre, F->mod, F->mod2);
        DFT2_NODE_LAZY_4_4(p0[k+1], p1[k+1], w, wpre, F->mod, F->mod2);
        DFT2_NODE_LAZY_4_4(p0[k+2], p1[k+2], w, wpre, F->mod, F->mod2);
        DFT2_NODE_LAZY_4_4(p0[k+3], p1[k+3], w, wpre, F->mod, F->mod2);
    }

    dft_node_lazy_4_4(p, depth-1, 2*node, F);
    n_fft_tft_draft(p1, len_rec, olen_rec, depth-1, 2*node+1, F);
}












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
void tft_node_lazy_4_4(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
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
        tft_node_lazy_4_4(p0, depth-2, 4*node, F);
        tft_node_lazy_4_4(p1, depth-2, 4*node+1, F);
        tft_node_lazy_4_4(p2, depth-2, 4*node+2, F);
        tft_node_lazy_4_4(p3, depth-2, 4*node+3, F);
    }
}

/** 2**depth-point DFT
 * Same specification as n_fft_tft, except for:
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 *             requirement (not checked): depth <= F.depth
 * * lazy_2_4: in [0..2n) / out [0..4n) / max < 4n
 *             requirement (not checked): 3 <= depth <= F.depth
 */
void tft_lazy_2_4(nn_ptr p, ulong depth, n_fft_args_t F)
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
        tft_lazy_2_4(p0, depth-2, F);
        tft_node_lazy_4_4(p1, depth-2, 1, F);
        tft_node_lazy_4_4(p2, depth-2, 2, F);
        tft_node_lazy_4_4(p3, depth-2, 3, F);
    }
}

void tft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F)
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
        tft_lazy_2_4(p0, depth-2, F);
        tft_node_lazy_4_4(p1, depth-2, 1, F);
        tft_node_lazy_4_4(p2, depth-2, 2, F);
        tft_node_lazy_4_4(p3, depth-2, 3, F);
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

void n_fft_tft(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth > 0)
    {
        n_fft_args_t Fargs;
        n_fft_set_args(Fargs, F->mod, F->tab_w);
        tft_lazy_1_4(p, depth, Fargs);
        for (ulong k = 0; k < (UWORD(1) << depth); k++)
        {
            if (p[k] >= Fargs->mod2)
                p[k] -= Fargs->mod2;
            if (p[k] >= Fargs->mod)
                p[k] -= Fargs->mod;
        }
    }
}

/*---------------*/
/* some comments */
/*---------------*/

/** TODO add any useful comment
 */

/** Base cases:
 * - having macros for "small" lengths (up to 16 or 32 at least) improves performance
 * - removing the base cases depth==3 in internal functions where this case is
 *   not really used (eg dft_node_lazy_4_4) does not make a difference
 */

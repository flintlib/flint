/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft.h"
#include "n_fft/impl.h"
#include "impl_macros_tft.h"
#include "ulong_extras.h"

/* division by x**d - c, lazy with precomputation */
/* coeff bounds: in [0, 4*n) | out [0, 4*n) */
/* TODO see if can be faster by using ideas from try_sparse */
FLINT_FORCE_INLINE
void _nmod_poly_divrem_circulant_lazy_4_4(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2)
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



/*-----------------------*/
/*  auxiliary functions  */
/*-----------------------*/

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

/* here ilen == len, len/2 < olen < len, olen multiple of 4 */
/* TODO TBD minimum len == 16, in which case olen == 12 */
void tft_node_lazy_4_4(nn_ptr p, ulong olen, ulong depth, ulong node, n_fft_args_t F)
{
    /* flint_printf("olen = %wu, len = %wu, node = %wu\n", olen, 1L<<depth, node); */
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
        ulong new_depth = n_clog2_gt2(olen);
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

/* ilen is arbitrary, > 0, multiple of 4 */
/* olen is arbitrary, > 0, multiple of 4 */
/* node == 0 */
/* p has space for max(ilen, len) coeffs where len is next power of 2 of olen */
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F)
{
    const ulong depth = n_clog2_gt2(olen);
    const ulong ilen_depth = n_clog2_gt2(ilen);
    if (ilen_depth < depth)  /* TODO special case for depth - 1? */
    {
        const ulong ilen2 = UWORD(1) << ilen_depth;
        const ulong nb_subcalls = UWORD(1) << (depth - ilen_depth);
        /* flint_printf("here %wu, %wu, %wu, %wu\n", olen, ilen_depth, depth, ilen2); */
        for (ulong k = ilen; k < ilen2; k++)
            p[k] = 0;
        for (ulong i = 1; i < nb_subcalls; i++)  /* TODO stop around olen */
        {
            for (ulong k = 0; k < ilen; k++)
                p[i * ilen2 + k] = p[k];
            for (ulong k = ilen; k < ilen2; k++)
                p[i * ilen2 + k] = 0;
        }
        ulong i;
        for (i = 0; olen >= ilen2 && i < nb_subcalls - 1; i++)
        {
            /* flint_printf("starting %wu, %wu, %wu, %wu --> %wu\n", olen, ilen_depth, depth, ilen2, i); */
            dft_node_lazy_4_4(p, ilen_depth, i, F);
            /* flint_printf("end of %wu, %wu, %wu, %wu --> %wu\n", olen, ilen_depth, depth, ilen2, i); */
            p = p + ilen2;
            olen -= ilen2;
        }
        /* flint_printf("here %wu, %wu, %wu, %wu --> %wu\n", olen, ilen_depth, depth, ilen2, depth - ilen_depth); */
        if (olen > 0)
            tft_node_lazy_4_4(p, olen, ilen_depth, i, F);
        /* flint_printf("alright %wu, %wu, %wu, %wu\n", olen, ilen_depth, depth, ilen2); */
    }
    else if (ilen_depth == depth)
    {
        const ulong len = UWORD(1) << ilen_depth;
        for (ulong k = ilen; k < len; k++)
            p[k] = 0;
        tft_node_lazy_4_4(p, olen, depth, 0, F);
    }
    else  /* ilen_depth > depth */
    {
        const ulong len = UWORD(1) << depth;
        _nmod_poly_divrem_circulant1(p, ilen, len, F->mod);
        tft_lazy_1_4(p, len, olen, F);
    }
}










/*-------------------*/
/*  main interfaces  */
/*-------------------*/

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

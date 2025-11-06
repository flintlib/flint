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
#include "n_fft/impl_macros_dft.h"
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

/* ilen is >= 8, multiple of 4 */
/* olen is >= 4, multiple of 4 */
/* exception: ilen==4 is accepted when olen==4 */
/* node == 0 */
/* p has space for max(ilen, len) coeffs where len is next power of 2 of olen */
/** NOTE: some calls could be lazier, but gain should be quite modest:
 *     [comment**]  DFT2_NODE_LAZY_4_4 -> LAZY_1_4
 *     [comment*]   dft_node_lazy_4_4 -> 2_4
 */
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F)
{
    /* base case */
    if (olen == 4 && ilen == 4)
    {
        DFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
        return;
    }

    const ulong odepth = n_clog2_gt2(olen);
    const ulong idepth = n_clog2_gt2(ilen);

    /* ilen >> olen: reduce mod x**len - 1 and call again */
    if (idepth > odepth)
    {
        const ulong len = UWORD(1) << odepth;
        _nmod_poly_divrem_circulant1(p, ilen, len, F->mod);
        tft_lazy_1_4(p, len, olen, F);
    }

    /* ilen ~ olen: recurse into dft + tft at half length */
    else if (idepth == odepth)
    {
        /* flint_printf("===equal=== %wu\n", odepth); */
        const ulong len = UWORD(1) << odepth;
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        ulong k;
        for (k = 0; k < ilen - len/2; k++)  /* FIXME try unrolling? */
            DFT2_LAZY_1_2(p0[k], p1[k], F->mod);
        for (; k < len/2; k++)
            p1[k] = p0[k];
        dft_lazy_2_4(p0, odepth - 1, F);
        tft_node_lazy_4_4(p1, olen - len/2, odepth - 1, 1, F);
    }

    /* if ilen ~ olen/2, do special butterfly and recurse */
    /* this is not necessary: it is covered by the next (convoluted) "if"; but */
    /* since this case is frequent, let's write it in a simple way */
    /* else if (idepth == odepth - 1) */
    /* { */

    /* } */

    else  /* idepth < odepth - 1 */
    {
        const ulong ilen2 = UWORD(1) << idepth;
        ulong i, k;
        ulong val;
        nn_ptr pp = p + ilen2;
        const ulong len_rec = ilen2 / 2;
        /* flint_printf("here %wu, %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2, len_rec); */

        /* first DFT round with node == 0 */
        /* save values and perform butterflies */
        for (k = 0; k < ilen - len_rec; k++)  /* FIXME try unrolling? */
        {
            pp[k] = p[k];
            pp[len_rec + k] = p[len_rec + k];
            DFT2_LAZY_1_2(p[k], p[len_rec + k], F->mod);
        }
        for (; k < len_rec; k++)
        {
            /* butterfly with input p[len_rec + k]: p[k] unchanged, p[len_rec + k] = p[k] */
            val = p[k];
            p[len_rec + k] = val;
            pp[k] = val;
        }
        /* flint_printf("before 2_4 %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */
        dft_lazy_2_4(p, idepth - 1, F);
        /* flint_printf("before 4_4 %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */
        dft_node_lazy_4_4(p+len_rec, idepth - 1, 1, F);  /* [comment*] */
        /* flint_printf("ok %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */

        /* update */
        p += ilen2;
        pp += ilen2;
        olen -= ilen2;
        /* flint_printf("here %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */

        /* next DFT rounds */
        for (i = 1; olen > ilen2; i++)  /* remove nb_subcalls? */
        {
            /* flint_printf("DFT rounds %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */
            /* save values and perform butterflies */
            for (k = 0; k < ilen - len_rec; k++)  /* FIXME try unrolling? */
            {
                pp[k] = p[k];
                pp[len_rec + k] = p[len_rec + k];
                DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment**] */
                                   F->tab_w[2*i], F->tab_w[2*i+1],
                                   F->mod, F->mod2);
            }
            for (; k < len_rec; k++)
            {
                /* butterfly with input p[len_rec + k]: p[k] unchanged, p[len_rec + k] = p[k] */
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

        /* flint_printf("final TFT round: %wu, %wu, %wu, %wu\n", olen, idepth, odepth, ilen2); */
        /* perform butterflies */
        for (k = 0; k < ilen - len_rec; k++)  /* FIXME try unrolling? */
        {
            DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment**] */
                               F->tab_w[2*i], F->tab_w[2*i+1],
                               F->mod, F->mod2);
        }
        for (; k < len_rec; k++)
            p[len_rec + k] = p[k];
        dft_node_lazy_4_4(p, idepth - 1, 2*i, F);
        /* flint_printf("arriving just before TFT\n"); */
        olen -= len_rec;
        if (olen > 0)
            tft_node_lazy_4_4(p+len_rec, olen, idepth - 1, 2*i+1, F);
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

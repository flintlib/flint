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
#include "nmod_poly.h"  // TODO probably not needed in the end
#include "nmod.h"  // TODO remove when poly_rem removed
#include "nmod_vec.h"

/* TODO temporary until something of this kind gets merged in */
/* (unless this version becomes too n_fft-specific and it will stay here) */
/* division by x**n - c, lazy with precomputation */
/* constraint: max(A) + 2*modn <= 2**FLINT_BITS */
/* coeff bounds: in [0, max(A)] | out [0, max(A) + 2*modn) */
/* TODO make faster by using ideas from try_sparse */
static inline void
_nmod_poly_divrem_xnmc_precomp_lazy(nn_ptr p, slong len, ulong n, ulong c, ulong c_precomp, ulong modn)
{
    /* assumes len >= n */
    slong i;
    ulong j, r, val, p_hi, p_lo;

    r = len % n;
    i = len - r - n;  /* multiple of n, >= 0 by assumption */

    for (j = 0; j < r; j++)
    {
        /* computes either val = (c*val mod n) or val = (c*val mod n) + n */
        val = p[i+n+j];
        umul_ppmm(p_hi, p_lo, c_precomp, val);
        val = c * val - p_hi * modn;
        /* lazy addition, yields RQ[i+j] in [0..k+2n), where max(RQ) <= k */
        p[i+j] = val + p[i+j];
    }

    i -= n;
    while (i >= 0)
    {
        for (j = 0; j < n; j++)
        {
            /* computes either val = (c*val mod n) or val = (c*val mod n) + n */
            val = p[i+n+j];
            umul_ppmm(p_hi, p_lo, c_precomp, val);
            val = c * val - p_hi * modn;
            /* lazy addition, yields RQ[i+j] in [0..k+2n), where max(RQ) <= k */
            p[i+j] = val + p[i+j];
        }
        i -= n;
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
void tft_node_lazy_4_4_v1(nn_ptr p, ulong ilen, ulong olen, ulong depth, ulong node, n_fft_args_t F)
{
    /* current full length */
    const ulong len = UWORD(1) << depth;

    /* if olen <= len/2, go down the tree */
    if (2*olen <= len)
    {
        /* flint_printf("HERE! 2*olen <= len, %wu, %wu\n", olen, len); */
        tft_node_lazy_4_4_v1(p, ilen, olen, depth-1, 2*node, F);
        return;
    }

    /* -> in what follows, len/2 < olen <= len */

    if (olen == 1)  /* also means len == 1: single point evaluation */
    {
        /* FIXME as such, this requires max(p) + 2n - 1 < 2**FLINT_BITS */
        /* FIXME efficiency: this evaluation should be done faster for point 1 or -1 */
        p[0] = _nmod_poly_evaluate_nmod_precomp_lazy(p, ilen, F->tab_w[2*node], F->tab_w[2*node+1], F->mod);
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

    /* if ilen <= len/2, directly do recursive calls */
    /* (this is a simplification of the general case below, and also */
    /* allows to do both calls with ilen again instead of len/2) */
    /* TODO think about changing this into a loop beforehand,
     * which would ensures we have len/2 < ilen <= len,
     * and then apply one butterfly which implies ilen == len?! */
    if (2*ilen <= len)
    {
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        for (ulong k = 0; k < ilen; k++)
            p1[k] = p0[k];
        tft_node_lazy_4_4_v1(p0, ilen,        len/2, depth-1, 2*node  , F);
        tft_node_lazy_4_4_v1(p1, ilen, olen - len/2, depth-1, 2*node+1, F);
        return;
    }

    /* if ilen > len, reduce f mod x**len - root, and continue */
    if (ilen > len)
    {
        /* once new function finalized, replace by:
         * _nmod_poly_divrem_xnmc_precomp_lazy(p, p, ilen, len, F->tab_w[2*node], F->tab_w[2*node+1], F->mod);
         * careful with lazy!! what guarantee on input? */
        nn_ptr R = flint_malloc(len * sizeof(ulong));
        nn_ptr B = flint_calloc(len+1, sizeof(ulong));  // x**n - root
        B[len] = 1;
        if (node % 2 == 0)
            B[0] = F->mod - F->tab_w[node];
        else  /* node % 2 == 1 */
            B[0] = F->mod - F->tab_w[node-1];
        nmod_t mod;
        nmod_init(&mod, F->mod);
        _nmod_poly_rem(R, p, ilen, B, len+1, mod);
        for (ulong k = 0; k < len; k++)
            p[k] = R[k];
        flint_free(R);
        flint_free(B);
        ilen = len;
    }

    /* TODO base cases to fix !! */
    if (len == olen)
    {
        /* flint_printf("base case len == olen, %wu, %wu, %wu\n", ilen, olen, len); */
        _nmod_vec_zero(p + ilen, len - ilen);
        dft_node_lazy_4_4(p, depth, node, F);
        return;
    }

    /* now len/2 < ilen <= len, do butterfly and then recursive calls */
    // in: [0..4n), out: [0..4n)
    const nn_ptr p0 = p;
    const nn_ptr p1 = p + len/2;
    const ulong w = F->tab_w[2*node];
    const ulong wpre = F->tab_w[2*node+1];
    /* TODO do this unroll 4 once base case ok */
    /* for (ulong k = 0; k < len/2; k+=4) */
    /* { */
    /*     DFT2_NODE_LAZY_4_4(p0[k+0], p1[k+0], w, wpre, F->mod, F->mod2); */
    /*     DFT2_NODE_LAZY_4_4(p0[k+1], p1[k+1], w, wpre, F->mod, F->mod2); */
    /*     DFT2_NODE_LAZY_4_4(p0[k+2], p1[k+2], w, wpre, F->mod, F->mod2); */
    /*     DFT2_NODE_LAZY_4_4(p0[k+3], p1[k+3], w, wpre, F->mod, F->mod2); */
    /* } */
    for (ulong k = 0; k < ilen - len/2; k++)
    {
        DFT2_NODE_LAZY_4_4(p0[k+0], p1[k+0], w, wpre, F->mod, F->mod2);
        // TODO for the moment we ensure p1 is fully reduced
        if (p1[k] >= F->mod2)
            p1[k] -= F->mod2;
        if (p1[k] >= F->mod)
            p1[k] -= F->mod;
    }
    for (ulong k = ilen - len/2; k < len/2; k++)
    {
        p1[k] = p0[k];
        /* // TODO for the moment we ensure p1 is fully reduced */
        if (p1[k] >= F->mod2)
            p1[k] -= F->mod2;
        if (p1[k] >= F->mod)
            p1[k] -= F->mod;
    }

    dft_node_lazy_4_4(p, depth-1, 2*node, F);
    tft_node_lazy_4_4_v1(p1, len/2, olen - len/2, depth-1, 2*node+1, F);
}

/* here ilen == len, len/2 < olen < len, olen multiple of 4 */
/* minimum len == 16, in which case olen == 12 */
void tft_node_lazy_4_4_v2_olen(nn_ptr p, ulong olen, ulong depth, ulong node, n_fft_args_t F)
{
    /* flint_printf("olen = %wu, len = %wu, node = %wu\n", olen, 1L<<depth, node); */
    if (depth == 2)
    {
        TFT4_2_NODE_LAZY_4_4(p[0], p[1], p[2], p[3],
                             F->tab_w[2*node], F->tab_w[2*node+1],
                             F->tab_w[4*node], F->tab_w[4*node+1],
                             F->mod, F->mod2);
    }
    else if (depth == 3)
    {
        TFT8_4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                             node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        if (olen == 4)
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
        else  /* olen == 16 */
        {
            TFT16_12_NODE_LAZY_4_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                                   p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                                   node, F->mod, F->mod2, F->tab_w);
        }
    }
    else if (depth == 5)
    {
        if (olen == 4)
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

        /* first recursive call in full length len/2 */
        dft_node_lazy_4_4(p0, depth-1, 2*node, F);

        // TODO at some point, ensure the relevant p1 is fully reduced
        const ulong olen_rec = olen - len/2;
        /* find out depth of call */
        /* e.g.: olen_rec == len/4+4 -> use (depth-1, 2*node+1) */
        /* e.g.: olen_rec == len/4-4 -> use (depth-2, 2*(2*node+1)) */
        depth -= 1;
        node = 2 * node + 1;
        while (olen_rec <= (UWORD(1) << (depth-1)))  /* computation to improve */
        {
            /* flint_printf("%wu ---- %wu, %wu\n", len/2, depth, node); */
            depth -= 1;
            node = 2 * node;
        }
        /* flint_printf("exiting loop with: %wu, %wu\n", depth, node); */
        /* now 1<<(depth - 1) < olen_rec < 1<<depth, */
        /* and newnode == (2 * oldnode + 1) << nb of while loop iterations */

        /* reduce p1 mod x**(len_rec) - root */
        ulong len_rec = UWORD(1) << depth;
        if (len_rec < len/2)
        {
            /* FIXME explain why this suffices; see if can be integrated in divrem */
            for (ulong k = 0; k < len/2; k++)
            {
                if (p1[k] >= F->mod2)
                    p1[k] -= F->mod2;
            }
            if (node % 2 == 0)
                _nmod_poly_divrem_xnmc_precomp_lazy(p1, len/2, len_rec, F->tab_w[node], F->tab_w[node+1], F->mod);
            else  /* node % 2 == 1 */
                _nmod_poly_divrem_xnmc_precomp_lazy(p1, len/2, len_rec, F->tab_w[node-1], F->tab_w[node], F->mod);
        }

        if (olen_rec == (UWORD(1) << depth))
        {
            if (olen_rec == 4)
            {
                DFT4_NODE_LAZY_4_4(p[0], p[1], p[2], p[3],
                                   F->tab_w[2*node], F->tab_w[2*node+1],
                                   F->tab_w[4*node], F->tab_w[4*node+1],
                                   F->tab_w[4*node+2], F->tab_w[4*node+3],
                                   F->mod, F->mod2);
            }
            else
                dft_node_lazy_4_4(p1, depth, node, F);
        }
        else
            tft_node_lazy_4_4_v2_olen(p1, olen_rec, depth, node, F);
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

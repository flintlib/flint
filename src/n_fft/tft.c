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

/** Truncated variants of the DFT provide finer control and performance
 * depending on the input polynomial length (ilen) and output number of
 * evaluations (olen).
 *
 * The TFT functions below take as input an array
 *    p = [p_0,...,p_{ilen-1},xx,...,xx]
 * of length max(ilen, len), where len is the smallest power of 2 greater than
 * or equal to olen. Here each `xx` can be any value, these values fill the
 * array up to the required length and might be modified during the algorithm.
 */


/*-----------------------*/
/*  auxiliary functions  */
/*-----------------------*/

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

/** truncated Fourier transform, general ilen, node 0
 * * similar to dft_lazy_1_4, but
 * accepts p of arbitrary length and only computes the first olen evaluations
 * * same requirements on depth and node as in dft_lazy_1_4, where here
 * depth is understood as ceiling(log_2(olen))
 * * Input requirements:
 *        ilen is a multiple of 4, at least 8
 *        olen is a positive multiple of 4
 *           (exception: ilen==4 is accepted when olen==4)
 *        p has space for max(ilen, len) coefficients where len = 1<<depth, and
 *        if ilen < len then coefficients ilen..len-1 may be modified by the
 *        algorithm
 * * lazy_1_4: in [0..n) / out [0..4n) / max < 4n
 */
/** NOTE: some calls could be lazier, but the gain should be quite modest:
 *     [comment**]  DFT2_NODE_LAZY_4_4 -> LAZY_1_4
 *     [comment*]   dft_node_lazy_4_4 -> 2_4
 */
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F)
{
    /* some recursive calls arrive at this base case */
    if (olen == 4 && ilen == 4)
    {
        DFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
        return;
    }

    const ulong odepth = n_clog2_gt2(olen);
    const ulong idepth = n_clog2_gt2(ilen);

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
                TFT2_1_NODE_LAZY_4_4(tmp0, tmp1, F->tab_w[2], F->tab_w[3], F->mod, F->mod2);  /* [comment**] */
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
                DFT2_NODE_LAZY_4_4(tmp0, tmp1, F->tab_w[2], F->tab_w[3], F->mod, F->mod2);  /* [comment**] */
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
        dft_node_lazy_4_4(p+len_rec, idepth - 1, 1, F);  /* [comment*] */

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
                DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment**] */
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
                DFT2_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment**] */
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
                TFT2_1_NODE_LAZY_4_4(p[k], p[len_rec + k],  /* [comment**] */
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


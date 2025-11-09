/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"  /* for precomp_shoup */
#include "n_fft.h"
#include "n_fft/impl.h"
#include "impl_macros_dft.h"

/** Structure.
 * - The main interface is n_fft_idft, it solves the problem at node 0
 *   (interpolating at all roots of unity of order 2**depth), as documented in
 *   n_fft.h.
 * - The core function is `idft_node_lazy_1_2`, which goes up the subproduct
 *   tree towards an arbitrary node in this tree; it takes input values in
 *   [0..n) and return values in [0..2n), following the idea of lazy
 *   butterflies highlighted by David Harvey [Faster arithmetic for
 *   number-theoretic transforms, Journal of Symbolic Computation, Volume 60,
 *   2014, pp 113-119]. This function does not scale the output by the inverse
 *   of 2**depth.
 * - This core function costs more than a iDFT at node 0, at least for small or
 *   smallish lengths. So a specific function for node 0 is given
 *   (`idft_lazy_1_4`), targeting input values in [0..n) and return values in
 *   [0..4n). The main function `n_fft_idft` just calls `idft_lazy_1_4`, and
 *   then scales the output value by the inverse of 2**depth, also ensuring the
 *   output is in [0..n).
 */

/*---------------------------*/
/* IDFT: auxiliary functions */
/*---------------------------*/

/** 2**depth-point inverse DFT, general node
 * * In-place transform p = [p[i] for 0 <= i < len], where len == 2**depth,
 * into the list of coefficients q = [q[j] for 0 <= j < len] of the unique
 * polynomial q(x) of degree < len such that p[i] == q(w[i])  for 0 <= i < len
 * * Here we write w[k] for 0 <= k < len/2, defined as
 *            w[2*k]   == F->tab_w[2**depth * node + 2*k]
 *            w[2*k+1] == - F->tab_w[2**depth * node + 2*k];
 * these are the len roots of the polynomial x**len - F->tab_w[2*node]
 * * Requirements (not checked):
 *        3 <= depth
 *        (node+1) * 2**depth < 2**F.depth (length of F->tab_w)
 * * lazy_1_2: in [0..n) / out [0..2n) / max < 4n
 */
void idft_node_lazy_1_2(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    if (depth == 3)
    {
        IDFT8_NODE_LAZY_1_2(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                            node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        IDFT16_NODE_LAZY_1_2(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                             p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                             node, F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        IDFT32_NODE_LAZY_1_2(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                             p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                             p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                             p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
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
        idft_node_lazy_1_2(p0, depth-2, 4*node, F);
        idft_node_lazy_1_2(p1, depth-2, 4*node+1, F);
        idft_node_lazy_1_2(p2, depth-2, 4*node+2, F);
        idft_node_lazy_1_2(p3, depth-2, 4*node+3, F);

        const ulong w2 = F->tab_w[2*node];
        const ulong w2_pr = F->tab_w[2*node+1];
        const ulong w = F->tab_w[4*node];
        const ulong w_pr = F->tab_w[4*node+1];
        const ulong Iw = F->tab_w[4*node+2];
        const ulong Iw_pr = F->tab_w[4*node+3];

        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_NODE_LAZY_2_2(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2);
            IDFT4_NODE_LAZY_2_2(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2);
            IDFT4_NODE_LAZY_2_2(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2);
            IDFT4_NODE_LAZY_2_2(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2);
        }
    }
}

/** 2**depth-point inverse DFT
 * Same specification as n_fft_idft, except that the
 * output values are in [0..4n)
 */
void idft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 0)
        return;

    if (depth == 1)
    {
        DFT2_LAZY_1_2(p[0], p[1], F->mod);
    }
    else if (depth == 2)
    {
        IDFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3],
                           F->mod, F->mod2);
    }
    else
    if (depth == 3)
    {
        IDFT8_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                       F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 4)
    {
        IDFT16_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                        p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                        F->mod, F->mod2, F->tab_w);
    }
    else if (depth == 5)
    {
        IDFT32_LAZY_1_4(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                        p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15],
                        p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23],
                        p[24], p[25], p[26], p[27], p[28], p[29], p[30], p[31],
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
        idft_lazy_1_4(p0, depth-2, F);
        idft_node_lazy_1_2(p1, depth-2, 1, F);
        idft_node_lazy_1_2(p2, depth-2, 2, F);
        idft_node_lazy_1_2(p3, depth-2, 3, F);

        // 4-point butterflies
        // input p0 in [0,4n), p1,p2,p3 in [0,2n)
        // output p0,p1,p2,p3 in [0,4n)
        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_LAZY_4222_4(p0[k+0], p1[k+0], p2[k+0], p3[k+0], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
            IDFT4_LAZY_4222_4(p0[k+1], p1[k+1], p2[k+1], p3[k+1], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
            IDFT4_LAZY_4222_4(p0[k+2], p1[k+2], p2[k+2], p3[k+2], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
            IDFT4_LAZY_4222_4(p0[k+3], p1[k+3], p2[k+3], p3[k+3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2);
        }
    }
}

/*---------------------------*/
/* ITFT: auxiliary functions */
/*---------------------------*/

//def itft2(f, iolen, depth, node, tab_w, tab_iw, scale=1, verbose=False):
//    """ inverse truncated Fourier transform
//    Input:
//        :p: coefficients [p_0,...,p_{iolen-1},xx,...,xx];  iolen <= 2**depth, 
//            xx,...,xx are any values to reach length NextPowerOfTwo(iolen)
//             (these values may be overwritten)
//        :iolen: nonnegative integer, number of total evaluations in input
//               and number of sought coefficients in output
//        :depth: nonnegative integer, current depth in recursive tree
//        :node: nonnegative integer, current root in recursive tree
//        :tab_w: powers [1,w**(d/2),w**(d/4),w**(3d/4),...] of all principal d-th
//                root of unity w, in bit reversed order, for some d power of 2
//                such that 2**depth <= d
//        :tab_iw: same list for the inverse 1/w
//    Requirements:
//        iolen <=  2**depth
//        (node+1) * 2**depth <= 2**F->depth (the length of tab_w and tab_iw)
//    Output:
//        in-place modify p so that the first iolen output coefficients p[:iolen]
//        are those of the unique polynomial p(x) of degree < iolen such that the
//        evaluation p(w_i) is equal to the input p[i]
//        where w_i = ...   (something close to tab_w[node * 2**depth + i] and its opposite)
//    Algo:
//        uses a direct reduction-tree approach
//        --> we are currently at node x**d - tab_w[node], where d = 2**depth
//        --> subtrees have as roots the nodes x**(d/2) - tab_w[2*node] and x**(d/2) - tab_w[2*node+1]
//    """
//    # iolen == 0: nothing to do
//    if iolen == 0:
//        return
//
//    # depth == 0 --> whatever f[0] is, keep it
//    if depth == 0:
//        f[0] *= scale
//        return
//
//    # now depth >= 1, d >= 2
//    d = 1 << depth
//
//    # if iolen <= d/2, go down the tree
//    if 2*iolen <= d:
//        # find NextPower2(iolen)
//        depth_rec = ceil(log(iolen, 2))
//        d_rec = 1 << depth_rec
//        d_quo = 1 << (depth - depth_rec)  # == d / d_rec
//        itft2(f[:d_rec], iolen, depth_rec, d_quo * node, tab_w, tab_iw, scale, verbose)
//        for k in range(iolen):
//            f[k] = d_quo * f[k]
//        return
//
//    # from here on, 1 <= d/2 < iolen <= d
//    d_rec = 1 << (depth-1)
//    iolen_rec = iolen - d_rec
//
//    # 1st call
//    itft2(f[:d_rec], d_rec, depth-1, 2*node, tab_w, tab_iw, scale, verbose)
//
//    # 2nd call
//    itft2(f[d_rec:], iolen_rec, depth-1, 2*node+1, tab_w, tab_iw, scale, verbose)
//
//    # deduce low and high parts
//    for k in range(iolen_rec):
//        tmp = f[k]
//        f[k] = tmp + f[d_rec+k]
//        f[d_rec+k] = tab_iw[2*node] * (tmp - f[d_rec+k])
//    for k in range(iolen_rec, d_rec):
//        # here f[d_rec + k] == 0
//        f[d_rec+k] = tab_iw[2*node] * f[k]
//    # efficiency: in the above loops, could we avoid/limit the multiplications by
//    # tab_iw[2*node] depending on what happens next in the reduction?
//    # OR reduce first a copy of the low part, before gathering both parts? (seems feasible but requiring a buffer?)
//
//    reduce_mod_prod_xnma(iolen, tab_w, node, depth, f, d)
//
//    return

/* TODO determine and explain lazy_x_x */
/* TODO think about base cases to support */
/* FIXME depth is useless? assumes len/2 < iolen <= len */
/* assumes iolen > 0 */
/* FIXME for the moment, assumes iolen multiple of 4. Think about it (more problematic for interpolation than for evaluation) */
void itft_node_lazy_x_x(nn_ptr p, ulong iolen, ulong depth, ulong node, n_fft_ctx_t F)
{
    if (depth == 1)
    {
        IDFT2_NODE_LAZY_2_2(p[0], p[1], F->tab_iw[2], F->tab_iw[3], F->mod, 2*F->mod);
    }

    else if (depth == 2)
    {
        IDFT4_NODE_LAZY_2_2(p[0], p[1], p[2], p[3],
                            F->tab_iw[2*node+0], F->tab_iw[2*node+1],    
                            F->tab_iw[4*node+0], F->tab_iw[4*node+1],   
                            F->tab_iw[4*node+2], F->tab_iw[4*node+3],
                            F->mod, 2*F->mod);
    }

    /* now depth >= 3, iolen >= 4 multiple of 4 */
    else
    {
        const ulong len = UWORD(1) << depth;

        /* iolen == len : call idft */
        if (iolen == len)
        {
            n_fft_args_t Fargs;
            n_fft_set_args(Fargs, F->mod, F->tab_iw);
            idft_node_lazy_1_2(p, depth, node, Fargs);
            return;
        }

        /* iolen <= len/2 : go down the tree */
        /* FIXME move this at the end for 2nd rec call; and here assume depth is the right one */
        if (2*iolen <= len)
        {
            ulong new_depth = n_clog2_ge2(iolen);  /* FIXME careful with iolen==1 */
            node = node << (depth - new_depth);
            ulong pow2 = UWORD(1) << (depth - new_depth);
            depth = new_depth;
            itft_node_lazy_x_x(p, iolen, depth, node, F);
            // FIXME pow2 = 1 << (depth - new_depth), computed above,
            // with precomputation -> store it in F?
            ulong pow2_pr = n_mulmod_precomp_shoup(pow2, F->mod);
            for (ulong k = 0; k < iolen; k++)  /* FIXME could unroll if iolen multiple of 4 */
            {
                /* p[k] = pow2 * p[k]; */
                /* FIXME lazy how, 1_2 ? */
                N_MULMOD_PRECOMP_LAZY(p[k], pow2, p[k], pow2_pr, F->mod);
            }
            return;
        }

        /* from here on, 1 <= len/2 < iolen <= len */
        /* full idft */
        const nn_ptr p0 = p;
        const nn_ptr p1 = p + len/2;
        n_fft_args_t Fargs;
        n_fft_set_args(Fargs, F->mod, F->tab_iw);
        idft_node_lazy_1_2(p0, depth-1, 2*node, Fargs);
        itft_node_lazy_x_x(p1, iolen - len/2, depth-1, 2*node+1, F);

        /* butterflies */
        ulong k = 0;
        for ( ; k < iolen - len/2; k++)
        {
            /* FIXME might benefit from more tmp? measure and see */
            IDFT2_NODE_LAZY_2_2(p0[k], p1[k],
                                Fargs->tab_w[2*node], Fargs->tab_w[2*node+1],
                                Fargs->mod, Fargs->mod2);
        }
        for ( ; k < len/2; k++)
        {
            /* here virtually p1[k] == 0 --> p1[k] = tab_iw[2*node] * p0[k] */
            /* FIXME lazy how, x_2 ? */
            N_MULMOD_PRECOMP_LAZY(p1[k], F->tab_iw[2*node], p0[k], F->tab_iw[2*node+1], F->mod);
        }
        /* notes for possible better performance: in the above loops, could we avoid/limit the multiplications by */
        /* tab_iw[2*node] depending on what happens next in the reduction? */
        /* OR reduce first a copy of the low part, before gathering both parts? (seems feasible but requiring a buffer?) */

        n_fft_set_args(Fargs, F->mod, F->tab_w);
        _nmod_poly_rem_prod_root1_lazy_4_4(p, len, iolen, depth, node, Fargs);
        for (ulong k = 0; k < iolen; k++)  /* TODO understand why this is necessary */
        {
            if (p[k] >= Fargs->mod2)
                p[k] -= Fargs->mod2;
        }
    }
}


/*-------------------*/
/*  main interfaces  */
/*-------------------*/

void n_fft_dft_t(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth > 0)
    {
        n_fft_args_t Fargs;
        n_fft_set_args(Fargs, F->mod, F->tab_w);
        idft_lazy_1_4(p, depth, Fargs);
        for (ulong k = 0; k < (UWORD(1) << depth); k++)
        {
            if (p[k] >= Fargs->mod2)
                p[k] -= Fargs->mod2;
            if (p[k] >= Fargs->mod)
                p[k] -= Fargs->mod;
        }
    }
}

void n_fft_idft(nn_ptr p, ulong depth, n_fft_ctx_t F)
{
    if (depth > 0)
    {
        n_fft_args_t Fargs;
        n_fft_set_args(Fargs, F->mod, F->tab_iw);
        idft_lazy_1_4(p, depth, Fargs);

        const ulong inv2 = F->tab_inv2[2*depth-2];
        const ulong inv2_pr = F->tab_inv2[2*depth-1];
        for (ulong k = 0; k < (UWORD(1) << depth); k++)
            p[k] = n_mulmod_shoup(inv2, p[k], inv2_pr, F->mod);
    }
}

void n_fft_itft(nn_ptr p, ulong iolen, n_fft_ctx_t F)
{
    if (iolen > 1)
    {
        ulong depth = n_clog2_ge2(iolen);
        /* n_fft_args_t Fargs; */
        /* n_fft_set_args(Fargs, F->mod, F->tab_w); */
        itft_node_lazy_x_x(p, iolen, depth, 0, F);
        const ulong inv2 = F->tab_inv2[2*depth-2];
        const ulong inv2_pr = F->tab_inv2[2*depth-1];
        for (ulong k = 0; k < iolen; k++)
            p[k] = n_mulmod_shoup(inv2, p[k], inv2_pr, F->mod);
    }
}

/*---------------*/
/* some comments */
/*---------------*/

/** In n_fft_idft, there is apparently no gain from using the lazy mulmod_shoup
 * variant whose output is in [0..2n) (so one may as well use the non-lazy one
 * which ensures output < n)              
 */

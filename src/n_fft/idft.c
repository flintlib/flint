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

// TODO doc
// TODO make sure this is tested (code coverage: including for small depths)
// larger base cases?
void idft_node_lazy_1_2(nn_ptr p, ulong depth, ulong node, n_fft_args_t F)
{
    //if (depth == 1)
    //{
    //    ulong p_hi, p_lo, tmp;
    //    IDFT2_NODE_LAZY_2_2(p[0], p[1], F->tab_w[2*node], F->tab_w[2*node+1], F->mod, F->mod2, p_hi, p_lo, tmp);
    //}
    //else if (depth == 2)
    //{
    //    ulong p_hi, p_lo;
    //    IDFT4_NODE_LAZY_1_2(p[0], p[1], p[2], p[3],
    //                      F->tab_w[2*node], F->tab_w[2*node+1],
    //                      F->tab_w[4*node], F->tab_w[4*node+1],
    //                      F->tab_w[4*node+2], F->tab_w[4*node+3],
    //                      F->mod, F->mod2, p_hi, p_lo);
    //}
    //else
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
        ulong p_hi, p_lo;

        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_NODE_LAZY_2_2(p0[k+0], p1[k+0], p2[k+0], p3[k+0], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE_LAZY_2_2(p0[k+1], p1[k+1], p2[k+1], p3[k+1], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE_LAZY_2_2(p0[k+2], p1[k+2], p2[k+2], p3[k+2], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
            IDFT4_NODE_LAZY_2_2(p0[k+3], p1[k+3], p2[k+3], p3[k+3], w2, w2_pr, w, w_pr, Iw, Iw_pr, F->mod, F->mod2, p_hi, p_lo);
        }
    }
}

void idft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F)
{
    if (depth == 0)
        return;

    if (depth == 1)
    {
        ulong tmp;
        DFT2_LAZY_1_2(p[0], p[1], F->mod, tmp);
    }
    else if (depth == 2)
    {
        ulong p_hi, p_lo;
        IDFT4_LAZY_1_4(p[0], p[1], p[2], p[3], F->tab_w[2], F->tab_w[3],
                           F->mod, F->mod2, p_hi, p_lo);
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
        ulong p_hi, p_lo;
        for (ulong k = 0; k < len/4; k+=4)
        {
            IDFT4_LAZY_4222_4(p0[k+0], p1[k+0], p2[k+0], p3[k+0], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY_4222_4(p0[k+1], p1[k+1], p2[k+1], p3[k+1], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY_4222_4(p0[k+2], p1[k+2], p2[k+2], p3[k+2], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
            IDFT4_LAZY_4222_4(p0[k+3], p1[k+3], p2[k+3], p3[k+3], F->tab_w[2], F->tab_w[3], F->mod, F->mod2, p_hi, p_lo);
        }
    }
}

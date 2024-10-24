/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "n_fft.h"

void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong depth, ulong p)
{
    if (depth < 3)
        depth = 3;
    if (max_depth < depth)
        depth = max_depth;

    // fill basic attributes
    F->mod = p;
    F->mod2 = 2*p;
    //F->mod4 = 4*p;
    F->max_depth = max_depth;
    F->depth = 3;  // to be able to call fit_depth below

    // fill tab_w2
    ulong pr_quo, pr_rem, ww;
    ww = w;
    n_mulmod_precomp_shoup_quo_rem(&pr_quo, &pr_rem, ww, p);
    F->tab_w2[2*(max_depth-2)] = ww;
    F->tab_w2[2*(max_depth-2)+1] = pr_quo;
    for (slong k = max_depth-3; k >= 0; k--)
    {
        // ww <- ww**2 and its precomputed quotient
        n_mulmod_and_precomp_shoup(&ww, &pr_quo, ww, ww, pr_quo, pr_rem, pr_quo, p);
        pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, p);
        F->tab_w2[2*k] = ww;
        F->tab_w2[2*k+1] = pr_quo;
    }
    // at this stage, pr_quo and pr_rem are for k == 0 i.e. for I == tab_w2[0]

    // fill tab_w for depth 3
    ulong len = UWORD(1) << (depth-1);  // len >= 4
    F->tab_w = (nn_ptr) flint_malloc(2*len * sizeof(ulong));

    F->tab_w[0] = UWORD(1);
    F->tab_w[1] = n_mulmod_precomp_shoup(UWORD(1), p);
    F->tab_w[2] = F->tab_w2[0];
    F->tab_w[3] = F->tab_w2[1];
    F->tab_w[4] = F->tab_w2[2];
    F->tab_w[5] = F->tab_w2[3];
    n_mulmod_and_precomp_shoup(F->tab_w+6, F->tab_w+7, F->tab_w2[0], F->tab_w2[2], pr_quo, pr_rem, F->tab_w2[3], p);

    // complete tab_w up to specified depth
    n_fft_ctx_fit_depth(F, depth);
}

void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p)
{
    FLINT_ASSERT(p > 2 && flint_clz(p) >= 3);    // 2 < p < 2**61
    FLINT_ASSERT(flint_ctz(p - UWORD(1)) >= 3);  // p-1 divisible by 8

    // find the constant and exponent such that p == c * 2**max_depth + 1
    const ulong max_depth = flint_ctz(p - UWORD(1));
    const ulong c = (p - UWORD(1)) >> max_depth;

    // find primitive root w of order 2**max_depth
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, c, p);

    // fill all attributes and tables
    n_fft_ctx_init2_root(F, w, max_depth, depth, p);
}

void n_fft_ctx_clear(n_fft_ctx_t F)
{
    flint_free(F->tab_w);
}

void n_fft_ctx_fit_depth(n_fft_ctx_t F, ulong depth)
{
    if (F->max_depth < depth)
        depth = F->max_depth;

    if (depth > F->depth)
    {
        ulong len = UWORD(1) << (depth-1);  // len >= 8 (since depth >= 4)
        F->tab_w = flint_realloc(F->tab_w, 2*len * sizeof(ulong));

        // tab_w[2] is w**(L/8) * tab_w[0], where L = 2**max_depth,
        // tab_w[2*4,2*6] is w**(L/16) * tab_w[2*0,2*2],
        // tab_w[2*8,2*10,2*12,2*14] is w**(L/32) * tab_w[2*0,2*2,2*4,2*6], etc.
        // recall tab_w2[2*d] == w**(L / 2**(d+2))
        ulong d = F->depth - 1;
        ulong llen = UWORD(1) << (F->depth-1);
        ulong ww, pr_quo, pr_rem;
        for ( ; llen < len; llen <<= 1, d += 1)
        {
            ww = F->tab_w2[2*d];
            pr_quo = F->tab_w2[2*d+1];
            pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, F->mod);
            // for each k, tab_w[2*(k+llen)] <- ww * tab_w[2*k], and deduce precomputation
            for (ulong k = 0; k < llen; k++)
                n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*k, F->tab_w + 2*llen + 2*k+1,
                                            ww, F->tab_w[2*k],
                                            pr_quo, pr_rem, F->tab_w[2*k+1], F->mod);
        }
        F->depth = depth;
    }
}

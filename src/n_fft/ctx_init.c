/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_fft.h"

/** Given the precomputed quotient a_pr for modular multiplication by a mod n,
 *          a_pr == floor(a * 2**FLINT_BITS / n)
 * where we assume 0 < a < n and n does not divide a * 2**FLINT_BITS,
 * this returns the quotient for mulmod by -a mod n,
 *          floor( (n-a) * 2**FLINT_BITS / n)
 *          == 2**FLINT_BITS - ceil(a * 2**FLINT_BITS / n)
 *          == 2**FLINT_BITS - a_pr
 *
 * Note: the requirement "n does not divide a * 2**FLINT_BITS" follows
 * from the other requirement 0 < a < n as soon as n is odd; in n_fft.h
 * we will only use this for odd primes
 */
FLINT_FORCE_INLINE ulong n_mulmod_precomp_shoup_negate(ulong a_pr)
{
    return UWORD_MAX - a_pr;
}

void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong cofactor, ulong depth, ulong p)
{
    if (depth < 3)
        depth = 3;
    if (max_depth < depth)
        depth = max_depth;

    // fill basic attributes
    F->mod = p;
    F->max_depth = max_depth;
    F->cofactor = cofactor;
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

    // fill tab_inv2
    for (ulong k = 0; k < max_depth; k++)
    {
        F->tab_inv2[2*k] = p - (cofactor << (max_depth - k-1));
        F->tab_inv2[2*k+1] = n_mulmod_precomp_shoup(F->tab_inv2[2*k], p);
    }

    // fill tab_w and tab_iw for depth 3
    ulong len = UWORD(1) << (depth-1);  // len >= 4
    F->tab_w = (nn_ptr) flint_malloc(2*len * sizeof(ulong));
    F->tab_iw = (nn_ptr) flint_malloc(2*len * sizeof(ulong));

    // w**0 == iw**0 == 1
    F->tab_w[0] = UWORD(1);
    F->tab_w[1] = n_mulmod_precomp_shoup(UWORD(1), p);
    F->tab_iw[0] = UWORD(1);
    F->tab_iw[1] = F->tab_w[1];

    // w**(L/4) == I and iw**(L/4) == -I, L == 2**max_depth
    F->tab_w[2] = F->tab_w2[0];
    F->tab_w[3] = F->tab_w2[1];
    F->tab_iw[2] = p - F->tab_w2[0];
    F->tab_iw[3] = n_mulmod_precomp_shoup_negate(F->tab_w2[1]);

    // w**(L/8) == J and w**(3L/8) == I*J
    F->tab_w[4] = F->tab_w2[2];
    F->tab_w[5] = F->tab_w2[3];
    n_mulmod_and_precomp_shoup(F->tab_w+6, F->tab_w+7, F->tab_w2[0], F->tab_w2[2], pr_quo, pr_rem, F->tab_w2[3], p);

    // iw**(L/8) == -I*J and iw**(3L/8) == -J
    F->tab_iw[4] = p - F->tab_w[6];
    F->tab_iw[5] = n_mulmod_precomp_shoup_negate(F->tab_w[7]);
    F->tab_iw[6] = p - F->tab_w[4];
    F->tab_iw[7] = n_mulmod_precomp_shoup_negate(F->tab_w[5]);

    // complete tab_w up to specified depth
    n_fft_ctx_fit_depth(F, depth);
}

void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p)
{
    FLINT_ASSERT(p > 2 && flint_clz(p) >= 2);    // 2 < p < 2**(FLINT_BITS-2)
    FLINT_ASSERT(flint_ctz(p - UWORD(1)) >= 3);  // p-1 divisible by 8

    // find the constant and exponent such that p == c * 2**max_depth + 1
    const ulong max_depth = flint_ctz(p - UWORD(1));
    const ulong cofactor = (p - UWORD(1)) >> max_depth;

    // find primitive root w of order 2**max_depth
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, cofactor, p);

    // fill all attributes and tables
    n_fft_ctx_init2_root(F, w, max_depth, cofactor, depth, p);
}

void n_fft_ctx_clear(n_fft_ctx_t F)
{
    flint_free(F->tab_w);
    flint_free(F->tab_iw);
}

void n_fft_ctx_fit_depth(n_fft_ctx_t F, ulong depth)
{
    if (F->max_depth < depth)
        depth = F->max_depth;

    if (depth > F->depth)
    {
        ulong len = UWORD(1) << (depth-1);  // len >= 8 (since depth >= 4)
        F->tab_w = flint_realloc(F->tab_w, 2*len * sizeof(ulong));
        F->tab_iw = flint_realloc(F->tab_iw, 2*len * sizeof(ulong));

        // tab_w[2] is w**(L/8) * tab_w[0], where L = 2**max_depth,
        // tab_w[2*4,2*6] is w**(L/16) * tab_w[2*0,2*2],
        // tab_w[2*8,2*10,2*12,2*14] is w**(L/32) * tab_w[2*0,2*2,2*4,2*6], etc.
        // recall tab_w2[2*k] == w**(L / 2**(k+2))
        ulong d = F->depth - 1;
        ulong llen = UWORD(1) << (F->depth-1);
        ulong ww, pr_quo, pr_rem;
        for ( ; llen < len; llen <<= 1, d += 1)
        {
            ww = F->tab_w2[2*d];
            pr_quo = F->tab_w2[2*d+1];
            pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, F->mod);
            // for each k, tab_w[2*(k+llen)] <- ww * tab_w[2*k], and deduce precomputation
            for (ulong k = 0; k < llen; k+=4)
            {
                n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*(k+0), F->tab_w + 2*llen + 2*(k+0)+1,
                                           ww, F->tab_w[2*(k+0)],
                                           pr_quo, pr_rem, F->tab_w[2*(k+0)+1], F->mod);
                n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*(k+1), F->tab_w + 2*llen + 2*(k+1)+1,
                                           ww, F->tab_w[2*(k+1)],
                                           pr_quo, pr_rem, F->tab_w[2*(k+1)+1], F->mod);
                n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*(k+2), F->tab_w + 2*llen + 2*(k+2)+1,
                                           ww, F->tab_w[2*(k+2)],
                                           pr_quo, pr_rem, F->tab_w[2*(k+2)+1], F->mod);
                n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*(k+3), F->tab_w + 2*llen + 2*(k+3)+1,
                                           ww, F->tab_w[2*(k+3)],
                                           pr_quo, pr_rem, F->tab_w[2*(k+3)+1], F->mod);

                F->tab_iw[2*llen + 2*(llen-1-(k+0))] = F->mod - F->tab_w[2*llen + 2*(k+0)];
                F->tab_iw[2*llen + 2*(llen-1-(k+0)) + 1] = n_mulmod_precomp_shoup_negate(F->tab_w[2*llen + 2*(k+0)+1]);
                F->tab_iw[2*llen + 2*(llen-1-(k+1))] = F->mod - F->tab_w[2*llen + 2*(k+1)];
                F->tab_iw[2*llen + 2*(llen-1-(k+1)) + 1] = n_mulmod_precomp_shoup_negate(F->tab_w[2*llen + 2*(k+1)+1]);
                F->tab_iw[2*llen + 2*(llen-1-(k+2))] = F->mod - F->tab_w[2*llen + 2*(k+2)];
                F->tab_iw[2*llen + 2*(llen-1-(k+2)) + 1] = n_mulmod_precomp_shoup_negate(F->tab_w[2*llen + 2*(k+2)+1]);
                F->tab_iw[2*llen + 2*(llen-1-(k+3))] = F->mod - F->tab_w[2*llen + 2*(k+3)];
                F->tab_iw[2*llen + 2*(llen-1-(k+3)) + 1] = n_mulmod_precomp_shoup_negate(F->tab_w[2*llen + 2*(k+3)+1]);
            }
        }

        F->depth = depth;
    }
}

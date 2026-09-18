/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "n_fft.h"

// return bit reversal index of k for given nbits:
// e.g. br_index([0,1,2,3], 4) == [0, 8, 4, 12]
static inline ulong br_index(ulong k, ulong nbits)
{
    k = ((k >> 1) & 0x55555555) | ((k & 0x55555555) << 1);
    k = ((k >> 2) & 0x33333333) | ((k & 0x33333333) << 2);
    k = ((k >> 4) & 0x0F0F0F0F) | ((k & 0x0F0F0F0F) << 4);
    k = ((k >> 8) & 0x00FF00FF) | ((k & 0x00FF00FF) << 8);
    k = ( k >> 16             ) | ( k               << 16);
#if FLINT_BITS == 64
    k = ( k >> 32             ) | ( k               << 32);
#endif // FLINT_BITS == 64

    return k >> (FLINT_BITS - nbits);
}

int test_one(n_fft_ctx_t F, ulong max_depth, ulong depth, ulong p, flint_rand_t state)
{
    // if depth < 3, init is supposed to behave as if depth == 3
    depth = FLINT_MAX(3, depth);

    // check all basic attributes
    if (F->mod != p)
        return 1;

    if (F->max_depth != max_depth)
        return 2;

    if ((1 + (F->cofactor << max_depth)) != p)
        return 3;

    if (F->depth != depth)
        return 4;

    // retrieve primitive root and its inverse
    const ulong w = F->tab_w2[2*(max_depth-2)];
    const ulong iw = n_invmod(w, p);

    // check the primitive root
    if (n_powmod2(w, UWORD(1)<<max_depth, p) != UWORD(1)
            || n_powmod2(w, UWORD(1)<<(max_depth-1), p) != p-UWORD(1))
        return 5;

    // check all entries of tab_w2
    for (ulong k = 0; k < max_depth-1; k++)
    {
        ulong w2 = F->tab_w2[2*k];
        if (w2 != n_powmod2(w, UWORD(1)<<(max_depth-2-k), p))
            return 6;
        if (F->tab_w2[2*k+1] != n_mulmod_precomp_shoup(w2, p))
            return 7;
    }

    // check all entries of tab_inv2
    for (ulong k = 0; k < max_depth; k++)
    {
        ulong inv2 = F->tab_inv2[2*k];
        if (inv2 != n_invmod((UWORD(1)<<(k+1)), p))
            return 8;
        if (F->tab_inv2[2*k+1] != n_mulmod_precomp_shoup(inv2, p))
            return 9;
    }

    // check a few random entries of tab_w and tab_iw
    for (ulong j = 0; j < 1000; j++)
    {
        ulong k = n_randint(state, UWORD(1) << (F->depth - 1));
        ulong exp = br_index(k, F->max_depth - 1);

        ulong wk = F->tab_w[2*k];
        if (wk != n_powmod2(w, exp, p))
            return 10;
        if (F->tab_w[2*k+1] != n_mulmod_precomp_shoup(wk, p))
            return 11;

        ulong iwk = F->tab_iw[2*k];
        if (iwk != n_powmod2(iw, exp, p))
            return 12;
        if (F->tab_iw[2*k+1] != n_mulmod_precomp_shoup(iwk, p))
            return 13;
    }

    return 0;
}

TEST_FUNCTION_START(n_fft_ctx_init2, state)
{
    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p, max_depth;
        if (i % 20 != 0)
        {
            // take random prime in [17, 2**(FLINT_BITS-2))
#if FLINT_BITS == 64
            ulong bits = 5 + n_randint(state, 58);
#else
            ulong bits = 5 + n_randint(state, 25);
#endif
            p = n_randprime(state, bits, 1);
            max_depth = flint_ctz(p-1);

            // we need p such that 8 divides p-1
            while (max_depth < 3)
            {
                p = n_randprime(state, bits, 1);
                max_depth = flint_ctz(p-1);
            }
        }
        else
        {
            // the above will most often have max_depth 3 or 4
            // every now and then we want p with larger max_depth
#if FLINT_BITS == 64
            max_depth = 40 + n_randint(state, 10);
#else
            max_depth = 10 + n_randint(state, 10);
#endif
            p = 1 + (UWORD(1) << max_depth);
            while (! n_is_prime(p))
                p += (UWORD(1) << max_depth);
            max_depth = flint_ctz(p-1);
        }

        // take depth between 0 and min(12, max_depth)
        ulong depth = n_randint(state, FLINT_MIN(12, max_depth));

        // init
        n_fft_ctx_t F;
        n_fft_ctx_init2(F, depth, p);
        
        int res = test_one(F, max_depth, depth, p, state);

        if (res)
            TEST_FUNCTION_FAIL(
                    "prime = %wu\n"
                    "root of unity = %wu\n"
                    "max_depth = %wu\n"
                    "depth = %wu\n"
                    "error code = %wu\n",
                    p, F->tab_w2[2*(max_depth-2)], max_depth, depth, res);

        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}

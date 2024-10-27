/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_H
#define N_FFT_H

#include "flint.h"

#define N_FFT_CTX_DEFAULT_DEPTH 12

#ifdef __cplusplus
extern "C" {
#endif

/**
 * TODO[short term] augment precomputations with inverse roots
 * TODO[short term] add testing for general variants, not only node0
 * TODO[long term] large depth can lead to heavy memory usage
 *              --> provide precomputation-free functions
 * TODO[long term] on zen4 (likely on other cpus as well) ctx_init becomes
 * slower at some point, losing a factor 4 or more, probably due to caching;
 * what is annoying is that the depth where it becomes slower is significantly
 * smaller (~13-14) when tab_iw has been incorporated compared to without
 * tab_iw (it was depth ~20-21); see if this can be understood, and maybe play
 * with vectorization for those simple functions
 * TODO[later] provide forward function which reduces output to [0..n) ?
 *    unclear this is useful... to be decided later
 */

/** n_fft context:
 * parameters and tabulated powers of the primitive root of unity "w".
 **/

typedef struct
{
    ulong mod;                 // modulus, odd prime
    ulong max_depth;           // maximum supported depth (w has order 2**max_depth)
    ulong depth;               // depth supported by current precomputation
    nn_ptr tab_w;              // tabulated powers of w, see below
    nn_ptr tab_iw;             // tabulated powers of 1/w, see below
    ulong tab_w2[128];         // powers w**(2**k), see below
} n_fft_ctx_struct;
typedef n_fft_ctx_struct n_fft_ctx_t[1];

/** Requirements (not checked upon init):
 *     - mod is an odd prime < 2**(FLINT_BITS-2)
 *     - max_depth must be >= 3 (so, 8 must divide mod - 1)
 * Total memory cost of precomputations for arrays tab_{w,iw}:
 *     at most 2 * (128 + 2**depth) ulong's
 */

/** tab_w2:
 *    - length 128, with undefined entries at index 2*max_depth and beyond
 *    - contains powers w**d for d a power of 2, and corresponding
 *    precomputations for modular multiplication:
 *       -- for 0 <= k < max_depth-1, tab_w2[2*k] = w**(2**(max_depth-2-k))
 *          and tab_w2[2*k+1] = floor(tab_w2[2*k] * 2**FLINT_BITS / mod)
 *       -- for 2*max_depth <= k < 128, tab_w2[k] is undefined
 *
 * --> one can retrieve w as tab_w2[2 * (max_depth-2)]
 * --> the first elements are tab_w2 = [I, I_pr, J, J_pr, ...]
 * where I is a square root of -1 and J is a square root of I
 */

/** tab_w:
 *     - length 2**depth
 *     - contains 2**(depth-1) first powers of w in (max_depth-1)-bit reversed order,
 *     and corresponding precomputations for modular multiplication:
 *        -- for 0 <= k < 2**(depth-1), tab_w[2*k] = w**(br[k])
 *           and tab_w[2*k+1] = floor(tab_w[2*k] * 2**FLINT_BITS / mod)
 *  where br = [0, 2**(max_depth-2), 2**(max_depth-3), 3 * 2**(max_depth-3), ...]
 *  is the bit reversal permutation of length 2**(max_depth-1)
 *  (https://en.wikipedia.org/wiki/Bit-reversal_permutation)
 *
 * In particular the first elements are
 *      tab_w = [1, 1_pr, I, I_pr, J, J_pr, IJ, IJ_pr, ...]
 * where I is a square root of -1, J is a square root of I, and IJ = I*J. Note
 * that powers of w beyond 2**(max_depth-1), for example -1, -I, -J, etc. are
 * not stored.
 **/

/** tab_iw: same as tab_w but for the primitive root 1/w */

/** Note for init functions, when depth is provided:
 *   - if it is < 3, it is pretended that it is 3
 *   - it it is more than F->max_depth (the maximum possible with the given
 *   prime), it is reduced to F->max_depth
 * After calling init, precomputations support DFTs of length up to 2**depth
 */

// initialize with given root and given depth
void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong depth, ulong mod);

// find primitive root, initialize with given depth
void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p);

// same, with default depth
FLINT_FORCE_INLINE
void n_fft_ctx_init_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong p)
{ n_fft_ctx_init2_root(F, w, max_depth, N_FFT_CTX_DEFAULT_DEPTH, p); }

FLINT_FORCE_INLINE
void n_fft_ctx_init(n_fft_ctx_t F, ulong p)
{ n_fft_ctx_init2(F, N_FFT_CTX_DEFAULT_DEPTH, p); }

// grows F->depth and precomputations to support DFTs of depth up to depth
void n_fft_ctx_fit_depth(n_fft_ctx_t F, ulong depth);

void n_fft_ctx_clear(n_fft_ctx_t F);


typedef struct
{
    ulong mod;                 // modulus, odd prime
    ulong mod2;                // 2*mod  (storing helps for speed)
    //ulong mod4;                // 4*mod  (storing helps for speed)
    nn_srcptr tab_w;           // tabulated powers of w, see below
} n_fft_args_struct;
typedef n_fft_args_struct n_fft_args_t[1];

FLINT_FORCE_INLINE
void n_fft_set_args(n_fft_args_t F, ulong mod, nn_srcptr tab_w)
{
    F->mod = mod;
    F->mod2 = 2*mod;
    F->tab_w = tab_w;
}





/** dft:
 * transforms / inverse transforms / transposed transforms
 * at length a power of 2
 */
void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F);
void n_fft_idft(nn_ptr p, ulong depth, n_fft_ctx_t F);  // TODO
void n_fft_dft_t(nn_ptr p, ulong depth, n_fft_ctx_t F);  // TODO (idft on inverted roots)
void n_fft_idft_t(nn_ptr p, ulong depth, n_fft_ctx_t F);  // TODO (dft on inverted roots)





#ifdef __cplusplus
}
#endif

#endif /* N_FFT_H */

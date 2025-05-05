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
 * TODO[short term] for perf vs simplifying code:
 *   - bench intermediate functions to make sure there is nothing surprising
 *   - check if the dft32 and idft32 macros actually help or not (remove them if not)
 *   - check if having these p_hi,p_lo tmp in macro args is useful or they can be removed
 * TODO[short term] add proper testing for inverse / transposed variants
 *
 * TODO[long term] large depth can lead to heavy memory usage
 *              --> provide precomputation-free functions
 * TODO[long term] on zen4 (likely on other cpus as well) ctx_init becomes
 * slower at some point, losing a factor 4 or more, probably due to caching;
 * what is annoying is that the depth where it becomes slower is significantly
 * smaller (~13-14) when tab_iw has been incorporated compared to without
 * tab_iw (it was depth ~20-21); see if this can be understood, and maybe play
 * with vectorization for those simple functions
 */


/*-------------------------------------------------*/
/* STRUCTURES FOR FFT CONTEXT / FUNCTION ARGUMENTS */
/*-------------------------------------------------*/


/** n_fft context:
 *     - basic parameters
 *     - precomputed powers of the primitive root of unity and its inverse
 *     - precomputed inverses of 2**k
 *
 * Requirements (not checked upon init):
 *     - mod is an odd prime < 2**(FLINT_BITS-2)
 *     - max_depth must be >= 3 (so, 8 must divide mod - 1)
 * Total memory cost of precomputations for arrays tab_{w,iw,w2,inv2}:
 *     at most 2 * (2*FLINT_BITS + 2**depth) ulong's
 *
 * For more details about the content of tab_{w,iw,w2,inv2}, see comments below
 **/
typedef struct
{
    ulong mod;                    // modulus, odd prime
    ulong max_depth;              // maximum supported depth (w has order 2**max_depth)
    ulong cofactor;               // prime is 1 + cofactor * 2**max_depth
    ulong depth;                  // depth supported by current precomputation
    nn_ptr tab_w;                 // precomputed powers of w
    nn_ptr tab_iw;                // precomputed powers of 1/w
    ulong tab_w2[2*FLINT_BITS];   // precomputed powers w**(2**k)
    ulong tab_inv2[2*FLINT_BITS]; // precomputed inverses of 2**k
} n_fft_ctx_struct;
typedef n_fft_ctx_struct n_fft_ctx_t[1];


/** n_fft arguments:
 *     - modulus mod
 *     - its double 2*mod (storing helps for speed)
 *     - precomputed powers of w
 * To be used as an argument in FFT functions. In some parts, providing this
 * instead of the whole context increased performance. Also, this facilitate
 * using the same function with both tab_w and tab_iw (by forming an fft_args
 * with Fargs->tab_w = F->tab_iw.
 **/
typedef struct
{
    ulong mod;                 // modulus, odd prime
    ulong mod2;                // 2*mod
    nn_srcptr tab_w;           // tabulated powers of w, see below
} n_fft_args_struct;
typedef n_fft_args_struct n_fft_args_t[1];


/** tab_w2:
 *    - length 2*FLINT_BITS, with undefined entries at index 2*(max_depth-1) and beyond
 *    - contains powers w**d for d a power of 2, and corresponding
 *    precomputations for modular multiplication:
 *       -- for 0 <= k < max_depth-1, tab_w2[2*k] = w**(2**(max_depth-2-k))
 *          and tab_w2[2*k+1] = floor(tab_w2[2*k] * 2**FLINT_BITS / mod)
 *       -- for 2*(max_depth-1) <= k < 2*FLINT_BITS, tab_w2[k] is undefined
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

/** tab_inv2:
 *     - length 2*FLINT_BITS, with undefined entries at index 2*max_depth and beyond
 *     - contains the modular inverses of 2**k, and corresponding
 *    precomputations for modular multiplication:
 *       -- for 0 <= k < max_depth, tab_inv2[2*k] = the inverse of 2**(k+1)
 *       modulo mod, and tab_inv2[2*k+1] = floor(tab_inv2[2*k] * 2**FLINT_BITS / mod)
 *       -- for 2*max_depth <= k < 2*FLINT_BITS, tab_inv2[k] is undefined
 * 
 * Recall F->mod == 1 + cofactor * 2**max_depth, so
 *          1 == F->mod - cofactor * 2**(max_depth - k) * 2**k
 * --> the inverse of 2**k in (0, F->mod) is
 *          F->mod - cofactor * 2**(max_depth - k),
 * we do not really need to store it, but we want the precomputations as well
 */


/*------------------------------------------*/
/* PRECOMPUTATIONS / CONTEXT INITIALIZATION */
/*------------------------------------------*/

/** Note for init functions, when depth is provided:
 *   - if it is < 3, it is pretended that it is 3
 *   - it it is more than F->max_depth (the maximum possible with the given
 *   prime), it is reduced to F->max_depth
 * After calling init, precomputations support DFTs of length up to 2**depth
 */

/* initialize with given root and given depth */
void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong cofactor, ulong depth, ulong mod);

/* find primitive root, initialize with given depth */
void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p);

/* same, with default depth */
FLINT_FORCE_INLINE
void n_fft_ctx_init_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong cofactor, ulong p)
{ n_fft_ctx_init2_root(F, w, max_depth, cofactor, N_FFT_CTX_DEFAULT_DEPTH, p); }

FLINT_FORCE_INLINE
void n_fft_ctx_init(n_fft_ctx_t F, ulong p)
{ n_fft_ctx_init2(F, N_FFT_CTX_DEFAULT_DEPTH, p); }

/* grows F->depth and precomputations to support DFTs of depth up to depth */
void n_fft_ctx_fit_depth(n_fft_ctx_t F, ulong depth);

void n_fft_ctx_clear(n_fft_ctx_t F);

FLINT_FORCE_INLINE
void n_fft_set_args(n_fft_args_t F, ulong mod, nn_srcptr tab_w)
{
    F->mod = mod;
    F->mod2 = 2*mod;
    F->tab_w = tab_w;
}

/*-----------------------------*/
/* DFT / IDFT / DFT_t / IDFT_t */
/*-----------------------------*/

/** forward and inverse transforms, and their transposes:
 *    - length is a power of 2, len == 2**depth
 *    - requirement of all functions (not checked): depth <= F.depth
 *    - the comments below describe algorithms that modify the input array p in
 *    place: in these comments p stands for the input p, whereas q stands
 *    for the array p after running the algorithm
 *    - below in comments we write w[k] for 0 <= k < len/2, defined as
 *            w[2*k]   == F->tab_w[2*k]
 *            w[2*k+1] == - F->tab_w[2*k]
 *    - hence the list w[k] for 0 <= k < len gives the len roots of the
 *    polynomial x**len - 1, which are all powers of the chosen len-th
 *    primitive root of unity, with exponents listed in bit reversed order
 *    - the matrix of DFT of length len is the len x len matrix
 *             DFT_{w,len} = [ w[i]**j ]_{0 <= i, j < len}
 */

/** dft: discrete Fourier transform (q = DFT_{w,len} * p)
 * In-place transform p = [p[j] for 0 <= j < len], seen as a polynomial p(x) of
 * degree < len, into its evaluations
 *     q == [p(w[i])  for 0 <= i < len],
 * where p(w[i]) = sum(p[j] * w[i]**j for 0 <= j < len)
 */

/** idft: inverse discrete Fourier transform (q = DFT_{w,len}^{-1} * p)
 * In-place transform p = [p[i] for 0 <= i < len] into the list of coefficients
 * q = [q[j] for 0 <= j < len] of the unique polynomial q(x) of degree < len
 * such that p[i] == q(w[i])  for 0 <= i < len
 */

/** dft_t: transposed discrete Fourier transform (q = p * DFT_{w,len})
 * In-place transform p = [p[i] for 0 <= i < len] into the list of weighted
 * power sums
 *        q == [PowerSum(p, w**j) for 0 <= j < len]
 * where PowerSum(p, w**j) == sum(p[i] * w[i]**j for 0 <= i < len)
 */

/** idft_t: transposed inverse discrete Fourier transform (q = p * DFT_{w,len}^{-1})
 * In-place transform p = [p[j] for 0 <= j < len] into the coefficients q =
 * [q[i] for 0 <= i < len] which appear in the partial fraction decomposition
 *      p(x) = sum_{0 <= i < len} q[i] / (1 - w[i] * x) + O(x**len)
 * where p(x) is the power series p(x) = sum_{0 <= j < len} p[j] x**j  + O(x**len)
 */

void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F);

void n_fft_idft(nn_ptr p, ulong depth, n_fft_ctx_t F);

void n_fft_dft_t(nn_ptr p, ulong depth, n_fft_ctx_t F);

void n_fft_idft_t(nn_ptr p, ulong depth, n_fft_ctx_t F);

#ifdef __cplusplus
}
#endif

#endif /* N_FFT_H */

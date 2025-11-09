/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_FFT_IMPL_H
#define N_FFT_IMPL_H

#include "longlong.h"  /* provides flint_clz */
#include "n_fft.h"

/*---------*/
/* helpers */
/*---------*/

/** Shoup's modular multiplication with precomputation, lazy
 * (does not perform the excess correction step)
 *  --> computes either r or r+n and store it is res, where r = (a*b) % n
 *  --> a_pr is the precomputation for n, p_hi and p_lo are temporaries
 */
#define N_MULMOD_PRECOMP_LAZY(res, a, b, a_pr, n)             \
do {                                                          \
    ulong p_hi, p_lo;                                         \
    umul_ppmm(p_hi, p_lo, (a_pr), (b));                       \
    res = (a) * (b) - p_hi * (n);                             \
} while(0)

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

/** n_fft arguments:
 *     - modulus mod
 *     - its double 2*mod (storing helps for speed)
 *     - precomputed powers of w
 * To be used as an argument in internal FFT functions. In some parts,
 * providing this instead of the whole context increased performance. Also,
 * this facilitate using the same function with both tab_w and tab_iw (by
 * forming an fft_args with Fargs->tab_w = F->tab_iw), see e.g. the
 * implementation of n_fft_idft.
 **/
typedef struct
{
    ulong mod;                 // modulus, odd prime
    ulong mod2;                // 2*mod
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

/* special divrems */
void _nmod_poly_divrem_circulant_lazy_4_4(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2);
void _nmod_poly_divrem_circulant_lazy_4_4_v0(nn_ptr p, slong len, ulong d, ulong c, ulong c_precomp, ulong n, ulong n2);
void _nmod_poly_divrem_circulant1(nn_ptr p, slong len, ulong d, ulong n);
void _nmod_poly_divrem_circulant1_v1(nn_ptr p, slong len, ulong d, ulong n);
void _nmod_poly_rem_circulant1(nn_ptr p, slong len, ulong d, ulong n);
void _nmod_poly_divrem_circulant_lazy_4_2_t(nn_ptr p, ulong len, ulong d, ulong c, ulong c_precomp, ulong n);
void _nmod_poly_divrem_circulant1_t(nn_ptr p, ulong len, ulong d);
void _nmod_poly_rem_prod_root1_lazy_4_4(nn_ptr p, ulong len, ulong d,
                                        ulong depth, ulong node, n_fft_args_t F);
void _nmod_poly_rem_prod_root1_t_lazy_4_4(nn_ptr p, ulong len, ulong d,
                                          ulong depth, ulong node, n_fft_args_t F);

/* dft internals */
void dft_node_lazy_4_4(nn_ptr p, ulong depth, ulong node, n_fft_args_t F);
void dft_lazy_2_4(nn_ptr p, ulong depth, n_fft_args_t F);
void dft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F);

/* idft internals */
void idft_node_lazy_1_2(nn_ptr p, ulong depth, ulong node, n_fft_args_t F);
void idft_lazy_1_4(nn_ptr p, ulong depth, n_fft_args_t F);

/* tft internals */
void tft_node_lazy_4_4(nn_ptr p, ulong olen, ulong depth, ulong node, n_fft_args_t F);
void tft_lazy_1_4(nn_ptr p, ulong ilen, ulong olen, n_fft_args_t F);

/* itft internals */
void itft_node_lazy_x_x(nn_ptr p, ulong iolen, ulong depth, ulong node, n_fft_ctx_t F);
/* FIXME see how to handle args! needs tab_w, tab_iw, tab_w2 */

/* some functions to get the right parameters (FIXME work in progress) */

/* exponent of next power of 2 for x >= 2 */
FLINT_FORCE_INLINE
ulong n_clog2_ge2(ulong x)
{
    return FLINT_BITS - flint_clz(x - 1);
}

FLINT_FORCE_INLINE
ulong _next_multiple_of_4(ulong x)
{
    return (x + 3) & UWORD(0xFFFFFFFFFFFFFFFC);
}


/* FIXME doc: get parameters for tft functions and ensure F "long enough" */
/* in: olen > 0, ilen >= 0 */
/* out: stores actual _ilen and _olen parameters to be used in tft calls, */
/*      and returns the required length for the input array */
FLINT_FORCE_INLINE
ulong n_fft_tft_prepare(ulong * _ilen, ulong * _olen, ulong ilen, ulong olen, n_fft_ctx_t F)
{
    /* FIXME check: no function should do anything nontrivial when ilen==0 */
    if (ilen == 0)
    {
        *_ilen = 0;
        *_olen = olen;
        return olen;
    }

    /* FIXME support this kind of base case or go directly to 4? */
    if (ilen <= 2 && olen <= 2)
    {
        *_ilen = 2;
        *_olen = 2;
        return 2;
    }

    const ulong odepth = n_clog2_ge2(olen);
    const ulong len = UWORD(1) << odepth;

    *_ilen = _next_multiple_of_4(ilen);
    *_olen = _next_multiple_of_4(olen);

    n_fft_ctx_fit_depth(F, odepth);

    return FLINT_MAX(*_ilen, len);
}

FLINT_FORCE_INLINE
ulong n_fft_itft_prepare(ulong * _iolen, ulong iolen, n_fft_ctx_t F)
{
    /* FIXME check: no function should do anything nontrivial when ilen==0 */
    if (iolen == 0)
    {
        *_iolen = 0;
        return iolen;
    }

    /* FIXME support this kind of base case or go directly to 4? */
    if (iolen <= 2)
    {
        *_iolen = 2;
        return 2;
    }

    const ulong iodepth = n_clog2_ge2(iolen);
    const ulong len = UWORD(1) << iodepth;

    *_iolen = _next_multiple_of_4(iolen);

    n_fft_ctx_fit_depth(F, iodepth);

    return len;
}

#endif  /* N_FFT_IMPL_H */


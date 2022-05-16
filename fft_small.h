/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FFT_SMALL_H
#define FFT_SMALL_H

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t
#include "flint.h"
#include "mpn_extras.h"

#define LG_BLK_SZ 8
#define BLK_SZ 256
#define BLK_SHIFT 10


#ifdef __cplusplus
 extern "C" {
#endif

FLINT_INLINE ulong n_pow2(int k)
{
    return UWORD(1) << k;
}

FLINT_INLINE ulong n_min(ulong a, ulong b)
{
    return FLINT_MIN(a, b);
}

FLINT_INLINE ulong n_max(ulong a, ulong b)
{
    return FLINT_MAX(a, b);
}

FLINT_INLINE ulong n_cdiv(ulong a, ulong b)
{
    return (a + b - 1)/b;
}

FLINT_INLINE ulong n_round_up(ulong a, ulong b)
{
    return n_cdiv(a, b)*b;
}

FLINT_INLINE ulong n_saturate_bits(ulong a) /* needs a better name */
{
    a |= a >> 1;
    a |= a >> 2;
    a |= a >> 4;
    a |= a >> 8;
    a |= a >> 16;
    return a;
}


FLINT_DLL void* flint_aligned_alloc(ulong alignment, ulong size);

FLINT_DLL void flint_aligned_free(void* p);

typedef struct {
    double* data;
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
    ulong depth;
    ulong blk_sz;
    double* w2s;
    ulong wtab_depth;
} sd_fft_ctx_struct;

typedef sd_fft_ctx_struct sd_fft_ctx_t[1];

FLINT_INLINE ulong sd_fft_ctx_offset(const sd_fft_ctx_t Q, ulong I)
{
    return (I << LG_BLK_SZ) + 4*(I >> (BLK_SHIFT+2));
}

FLINT_INLINE ulong sd_fft_ctx_data_size(const sd_fft_ctx_t Q)
{
    return sd_fft_ctx_offset(Q, n_pow2(Q->depth - LG_BLK_SZ));
}

FLINT_INLINE void sd_fft_ctx_set_data(sd_fft_ctx_t Q, double* d)
{
    Q->data = d;
}

FLINT_INLINE double* sd_fft_ctx_release_data(sd_fft_ctx_t Q)
{
    double* d = Q->data;
    Q->data = NULL;
    return d;
}

FLINT_INLINE double* sd_fft_ctx_blk_index(const sd_fft_ctx_t Q, ulong I)
{
    return Q->data + sd_fft_ctx_offset(Q, I);
}

FLINT_INLINE void sd_fft_ctx_set_index(const sd_fft_ctx_t Q, ulong i, double x)
{
    sd_fft_ctx_blk_index(Q, i/BLK_SZ)[i%BLK_SZ] = x;
}

FLINT_INLINE double sd_fft_ctx_get_index(const sd_fft_ctx_t Q, ulong i)
{
    return sd_fft_ctx_blk_index(Q, i/BLK_SZ)[i%BLK_SZ];
}

FLINT_INLINE double sd_fft_ctx_get_fft_index(const sd_fft_ctx_t Q, ulong i)
{
    ulong j = i&(BLK_SZ-16);
    FLINT_ASSERT(BLK_SZ >= 32);
    j |= (i&3)<<2;
    j |= ((i>>2)&3);
    return sd_fft_ctx_blk_index(Q, i/BLK_SZ)[j];
}

FLINT_INLINE void sd_fft_ctx_clear(sd_fft_ctx_t Q)
{
    flint_aligned_free(Q->w2s);
}

/* sd_fft.c */
FLINT_DLL void sd_fft_base(const sd_fft_ctx_t Q, ulong I, ulong j);
FLINT_DLL void sd_fft_main_block(const sd_fft_ctx_t Q, ulong I, ulong S, ulong k, ulong j);
FLINT_DLL void sd_fft_main(const sd_fft_ctx_t Q, ulong I, ulong S, ulong k, ulong j);
FLINT_DLL void sd_fft_trunc_block(const sd_fft_ctx_t Q, ulong I, ulong S, ulong k, ulong j, ulong itrunc, ulong otrunc);
FLINT_DLL void sd_fft_trunc(const sd_fft_ctx_t Q, ulong I, ulong S, ulong k, ulong j, ulong itrunc, ulong otrunc);

/* sd_ifft.c */
void sd_ifft_trunc(
    const sd_fft_ctx_t Q,
    ulong I, // starting index
    ulong S, // stride
    ulong k, // transform length 2^(k + LG_BLK_SZ)
    ulong j,
    ulong z,   // actual trunc is z*BLK_SZ
    ulong n,   // actual trunc is n*BLK_SZ
    int f);

/* sd_fft_ctx.c */
FLINT_DLL void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp);
FLINT_DLL void sd_fft_ctx_fit_wtab(sd_fft_ctx_t Q, ulong k);
FLINT_DLL void sd_fft_ctx_set_depth(sd_fft_ctx_t Q, ulong l);

FLINT_INLINE void sd_fft_ctx_fft_trunc(const sd_fft_ctx_t Q, ulong itrunc, ulong otrunc)
{
    FLINT_ASSERT(itrunc % BLK_SZ == 0);
    FLINT_ASSERT(otrunc % BLK_SZ == 0);
    sd_fft_trunc(Q, 0, 1, Q->depth - LG_BLK_SZ, 0, itrunc/BLK_SZ, otrunc/BLK_SZ);
}

FLINT_INLINE void sd_fft_ctx_ifft_trunc(const sd_fft_ctx_t Q, ulong trunc)
{
    FLINT_ASSERT(trunc % BLK_SZ == 0);
    sd_ifft_trunc(Q, 0, 1, Q->depth - LG_BLK_SZ, 0, trunc/BLK_SZ, trunc/BLK_SZ, 0);
}

#ifdef __cplusplus
}
#endif

#endif


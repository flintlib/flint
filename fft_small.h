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
#include "machine_vectors.h"

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

/* 0 -> 0, 1 -> 1, [2,3] -> 3, [4,7] -> 7, [8,15] -> 15, ... */
FLINT_INLINE ulong n_next_pow2m1(ulong a)
{
    a |= a >> 1;
    a |= a >> 2;
    a |= a >> 4;
    a |= a >> 8;
    a |= a >> 16;
#if FLINT64
    a |= a >> 32;
#endif
    return a;
}

#if 0
FLINT_INLINE ulong n_clog2(ulong x) {
    if (x <= 2)
        return x == 2;

   ulong zeros = FLINT_BITS;
   count_leading_zeros(zeros, x - 1);
   return FLINT_BITS - zeros;
}
#endif


FLINT_INLINE ulong n_leading_zeros(ulong x) {
    return __builtin_clzll(x);
}

FLINT_INLINE ulong n_trailing_zeros(ulong x) {
    return __builtin_ctzll(x);
}

FLINT_INLINE ulong n_nbits(ulong x) {
    return 64 - n_leading_zeros(x);
}

FLINT_INLINE ulong n_clog2(ulong x) {
    return (x <= 2) ? (x == 2) : 64 - __builtin_clzll(x - 1);
}


FLINT_DLL void* flint_aligned_alloc(ulong alignment, ulong size);

FLINT_DLL void flint_aligned_free(void* p);

/*
    The twiddle factors are split across FLINT_BITS tables:

        [0] = {e(1)}                                original index 0
        [1] = {e(1/4)}                              original index 1
        [2] = {e(1/8), e(3/8)}                      original index 2,3
        [3] = {e(1/16), e(5/16), e(3/16), e(7/16)}  original index 4,5,6,7
        ...

    The unallocated ones start out as NULL, and once they are filled in, they
    never have to move. This simplifies thread safe enlargement but complicates
    random access into the original table. If j is the index into the original
    table, the new indicies are

        [j_bits][j_r]  where j_bits = nbits(j), j_r = j - 2^(j_bits-1)

    with the special case j_bits = j_r = 0 for j = 0.
    The first SD_FFT_CTX_INIT_DEPTH tables are stored consecutively to ease the
    lookup of small indices, which must currently be at least 4.
*/

#define SD_FFT_CTX_INIT_DEPTH 10

/* This context is the one expected to sit in a global position */
typedef struct {
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
    ulong blk_sz;
    volatile ulong w2tab_depth;
    double* w2tab[FLINT_BITS];
} sd_fft_ctx_struct;

typedef sd_fft_ctx_struct sd_fft_ctx_t[1];

/* The local context is expected to be copied and passed to the calculations. */
typedef struct {
    double p;
    double pinv;
    const double* w2tab[50];
} sd_fft_lctx_struct;

typedef sd_fft_lctx_struct sd_fft_lctx_t[1];


FLINT_INLINE ulong sd_fft_ctx_blk_offset(ulong I)
{
    return (I << LG_BLK_SZ) + 4*(I >> (BLK_SHIFT+2));
}

FLINT_INLINE ulong sd_fft_ctx_data_size(ulong depth)
{
    return sd_fft_ctx_blk_offset(n_pow2(depth - LG_BLK_SZ));
}

FLINT_INLINE double* sd_fft_ctx_blk_index(double* d, ulong I)
{
    return d + sd_fft_ctx_blk_offset(I);
}

FLINT_INLINE void sd_fft_ctx_set_index(double* d, ulong i, double x)
{
    sd_fft_ctx_blk_index(d, i/BLK_SZ)[i%BLK_SZ] = x;
}

FLINT_INLINE double sd_fft_ctx_get_index(double* d, ulong i)
{
    return sd_fft_ctx_blk_index(d, i/BLK_SZ)[i%BLK_SZ];
}

/* slightly-worse-than-bit-reversed order */
FLINT_INLINE double sd_fft_ctx_get_fft_index(double* d, ulong i)
{
    ulong j = i&(BLK_SZ-16);
    FLINT_ASSERT(BLK_SZ >= 16);
    j |= (i&3)<<2;
    j |= ((i>>2)&3);
    return sd_fft_ctx_blk_index(d, i/BLK_SZ)[j];
}

/* sd_fft.c */
FLINT_DLL void sd_fft_trunc(const sd_fft_lctx_t Q, double* d, ulong I, ulong S, ulong k, ulong j, ulong itrunc, ulong otrunc);

/* sd_ifft.c */
FLINT_DLL void sd_ifft_trunc(const sd_fft_lctx_t Q, double* d, ulong I, ulong S, ulong k, ulong j, ulong z, ulong n, int f);

/* sd_fft_ctx.c */
FLINT_DLL void sd_fft_ctx_clear(sd_fft_ctx_t Q);
FLINT_DLL void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp);
FLINT_DLL void sd_fft_ctx_fit_depth(sd_fft_ctx_t Q, ulong k);

FLINT_INLINE void sd_fft_ctx_fft_trunc(sd_fft_ctx_t Q, double* d, ulong depth, ulong itrunc, ulong otrunc)
{
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    FLINT_ASSERT(itrunc % BLK_SZ == 0);
    FLINT_ASSERT(otrunc % BLK_SZ == 0);
    sd_fft_ctx_fit_depth(Q, depth);
    sd_fft_lctx_t QL;
    QL->p = Q->p;
    QL->pinv = Q->pinv;
    for (int i = 0; i < 50; i++)
        QL->w2tab[i] = Q->w2tab[i];
    sd_fft_trunc(QL, d, 0, 1, depth - LG_BLK_SZ, 0, itrunc/BLK_SZ, otrunc/BLK_SZ);
}

FLINT_INLINE void sd_fft_ctx_ifft_trunc(sd_fft_ctx_t Q, double* d, ulong depth, ulong trunc)
{
    FLINT_ASSERT(depth >= LG_BLK_SZ);
    FLINT_ASSERT(trunc % BLK_SZ == 0);
    sd_fft_ctx_fit_depth(Q, depth);
    sd_fft_lctx_t QL;
    QL->p = Q->p;
    QL->pinv = Q->pinv;
    for (int i = 0; i < 50; i++)
        QL->w2tab[i] = Q->w2tab[i];
    sd_ifft_trunc(QL, d, 0, 1, depth - LG_BLK_SZ, 0, trunc/BLK_SZ, trunc/BLK_SZ, 0);
}

FLINT_INLINE double sd_fft_lctx_w2s(const sd_fft_lctx_t Q, ulong j)
{
    if (j == 0)
        return Q->w2tab[0][0];
    ulong n = n_nbits(j);
    return Q->w2tab[n][j-n_pow2(n-1)];
}


FLINT_INLINE double sd_fft_ctx_w2s(const sd_fft_ctx_t Q, ulong j)
{
    if (j == 0)
        return Q->w2tab[0][0];
    ulong n = n_nbits(j);
    return Q->w2tab[n][j-n_pow2(n-1)];
}

typedef struct {
    ulong prime;
    ulong coeff_len;
    ulong nprimes;
    ulong* data;
} crt_data_struct;

typedef crt_data_struct crt_data_t[1];

typedef void (*to_ffts_func)(
        sd_fft_ctx_struct* Qffts, double* d, ulong dstride,
        const ulong* a_, ulong an_, ulong atrunc,
        const vec4d* two_pow);

typedef void (*from_ffts_func)(
        ulong* z, ulong zn, ulong zlen,
        sd_fft_ctx_struct* Qffts, double* d, ulong dstride,
        crt_data_struct* Qcrts,
        ulong bits);

typedef struct {
    ulong np;
    ulong bits;
    ulong bn_bound;
    to_ffts_func to_ffts;
    from_ffts_func from_ffts;
} profile_entry_struct;

typedef profile_entry_struct profile_entry_t[1];

#define MPN_CTX_NSLOTS 8
#define MAX_NPROFILES 40
#define VEC_SZ 4

typedef struct {
    vec4d* data;
    ulong length;
} vec4dptr_with_length;

typedef struct {
    sd_fft_ctx_struct ffts[MPN_CTX_NSLOTS];
    crt_data_struct crts[MPN_CTX_NSLOTS];
    vec4dptr_with_length two_powers[MPN_CTX_NSLOTS];
    profile_entry_struct profiles[MAX_NPROFILES];
    ulong profiles_size;
    void* buffer;
    ulong buffer_alloc;
} mpn_ctx_struct;

typedef mpn_ctx_struct mpn_ctx_t[1];

void mpn_ctx_init(mpn_ctx_t R, ulong p);
void mpn_ctx_clear(mpn_ctx_t R);
void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong* z, ulong* a, ulong an, ulong* b, ulong bn);

#ifdef __cplusplus
}
#endif

#endif


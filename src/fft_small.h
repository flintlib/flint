/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FFT_SMALL_H
#define FFT_SMALL_H

#include "longlong.h"
#include "machine_vectors.h"

#if FLINT_USES_PTHREAD
# include <pthread.h>
# include <stdatomic.h>
#endif

#define LG_BLK_SZ 8
#define BLK_SZ 256

#ifdef __cplusplus
extern "C" {
#endif

/* Check that a modulus n satisfies the assumptions for mulmod
   documented in machine_vectors.h */
int fft_small_mulmod_satisfies_bounds(ulong n);

FLINT_FORCE_INLINE ulong n_pow2(int k)
{
    return UWORD(1) << k;
}

FLINT_FORCE_INLINE ulong n_min(ulong a, ulong b)
{
    return FLINT_MIN(a, b);
}

FLINT_FORCE_INLINE ulong n_max(ulong a, ulong b)
{
    return FLINT_MAX(a, b);
}

FLINT_FORCE_INLINE ulong n_cdiv(ulong a, ulong b)
{
    /* not technically correct because the addition can overflow */
    return (a + b - 1)/b;
}

FLINT_FORCE_INLINE ulong n_round_up(ulong a, ulong b)
{
    return n_cdiv(a, b)*b;
}

FLINT_FORCE_INLINE ulong n_round_down(ulong a, ulong b)
{
    return a/b*b;
}

/* 0 -> 0, 1 -> 1, [2,3] -> 3, [4,7] -> 7, [8,15] -> 15, ... */
FLINT_FORCE_INLINE ulong n_next_pow2m1(ulong a)
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

FLINT_FORCE_INLINE ulong n_leading_zeros(ulong x) {
    return x == 0 ? FLINT_BITS : flint_clz(x);
}

FLINT_FORCE_INLINE ulong n_trailing_zeros(ulong x) {
    return x == 0 ? FLINT_BITS : flint_ctz(x);
}

/*
    nbits is a mess

    without assuming x != 0:
        on x86 we want 64 - LZCNT
        on arm we want 64 - CLZ

        the problem is gcc decided that __builtin_clz is undefined on zero
        input even though both instructions LZCNT and CLZ are defined

    assuming x != 0:
        on x86 we want BSR + 1
*/
FLINT_FORCE_INLINE ulong n_nbits(ulong x) {
    if (x == 0)
        return 0;
    return FLINT_BITS - flint_clz(x);
}

FLINT_FORCE_INLINE ulong n_nbits_nz(ulong x) {
    FLINT_ASSERT(x != 0);
    return (flint_clz(x)^(FLINT_BITS-1)) + 1;
}

FLINT_FORCE_INLINE ulong n_clog2(ulong x) {
    return (x <= 2) ? (x == 2) : FLINT_BITS - flint_clz(x - 1);
}

FLINT_FORCE_INLINE ulong n_flog2(ulong x) {
    return (x <= 2) ? (x == 2) : FLINT_BITS - flint_clz(x);
}

FLINT_FORCE_INLINE slong z_min(slong a, slong b) {return FLINT_MIN(a, b);}

FLINT_FORCE_INLINE slong z_max(slong a, slong b) {return FLINT_MAX(a, b);}

int fft_small_mulmod_satisfies_bounds(ulong n);

/*
    The twiddle factors are split across SD_FFT_CTX_W2TAB_SIZE tables:

        [0] = {e(1)}                                original index 0
        [1] = {e(1/4)}                              original index 1
        [2] = {e(1/8), e(3/8)}                      original index 2,3
        [3] = {e(1/16), e(5/16), e(3/16), e(7/16)}  original index 4,5,6,7
        ...

    The unallocated ones start out as NULL, and once they are filled in, they
    never have to move. This simplifies thread safe enlargement but complicates
    random access into the original table. If j is the index into the original
    table, the new indices are

        [j_bits][j_r]  where j_bits = nbits(j), j_r = j - 2^(j_bits-1)

    with the special case j_bits = j_r = 0 for j = 0.
    The first SD_FFT_CTX_W2TAB_INIT tables are stored consecutively to ease the
    lookup of small indices, which must currently be at least max(4, LG_BLK_SZ).
*/

/* for the fft look up of powers of w */
#define SET_J_BITS_AND_J_R(j_bits, j_r, j) \
do { \
    if (j == 0) \
    { \
        j_bits = 0; \
        j_r = 0; \
    } \
    else \
    { \
        j_bits = n_nbits_nz(j); \
        j_r = j - n_pow2(j_bits - 1); \
    } \
} while (0)

/* for the ifft look up of powers of w^-1: the remainder has to be flipped */
#define SET_J_BITS_AND_J_MR(j_bits, j_mr, j) \
do { \
    if (j == 0) \
    { \
        j_bits = 0; \
        j_mr = 0; \
    } \
    else \
    { \
        j_bits = n_nbits_nz(j); \
        j_mr = n_pow2(j_bits) - 1 - j; \
    } \
} while (0)


#define SD_FFT_CTX_W2TAB_INIT 12
#define SD_FFT_CTX_W2TAB_SIZE 40

/* This context is the one expected to sit in a global position */
typedef struct {
    double p;
    double pinv;
    nmod_t mod;
    ulong primitive_root;
#if FLINT_USES_PTHREAD
    _Atomic(unsigned int) w2tab_depth;
#else
    unsigned int w2tab_depth;
#endif
    double* w2tab[SD_FFT_CTX_W2TAB_SIZE];
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
} sd_fft_ctx_struct;

typedef sd_fft_ctx_struct sd_fft_ctx_t[1];

FLINT_FORCE_INLINE ulong sd_fft_ctx_blk_offset(ulong I)
{
    return I << LG_BLK_SZ;
}

FLINT_FORCE_INLINE ulong sd_fft_ctx_data_size(ulong L)
{
    return n_pow2(L);
}

FLINT_FORCE_INLINE double* sd_fft_ctx_blk_index(double* d, ulong I)
{
    return d + sd_fft_ctx_blk_offset(I);
}

/*
location of the bit-reversed eval:
    with out_data = fft(in_data) of length 2^L, then
    eval_poly(in_data, sd_fft_ctx_w(, i)) = out_data[sd_fft_ctx_trunc_index(L, i)]
*/
FLINT_FORCE_INLINE ulong sd_fft_ctx_trunc_index(ulong L, ulong i)
{
    /* 4x4 transposed blocks in basecases if depth >= 4 */
    if (L >= 4)
        i = (i&(-16)) | ((i>>2)&3) | ((i&3)<<2);
    return i;
}

/* sd_fft.c */
void sd_fft_trunc(sd_fft_ctx_t Q, double* d, ulong L, ulong itrunc, ulong otrunc);

/* sd_ifft.c */
void sd_ifft_trunc(sd_fft_ctx_t Q, double* d, ulong L, ulong trunc);

/* sd_fft_ctx.c */
void sd_fft_ctx_clear(sd_fft_ctx_t Q);
void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp);
void sd_fft_ctx_fit_depth_with_lock(sd_fft_ctx_t Q, ulong k);

FLINT_FORCE_INLINE void sd_fft_ctx_fit_depth(sd_fft_ctx_t Q, ulong depth)
{
#if FLINT_USES_PTHREAD
    ulong tdepth = (ulong)atomic_load_explicit(&Q->w2tab_depth, memory_order_relaxed);
#else
    ulong tdepth = (ulong)Q->w2tab_depth;
#endif
    if (FLINT_UNLIKELY(tdepth < depth))
        sd_fft_ctx_fit_depth_with_lock(Q, depth);
}

void sd_fft_ctx_point_mul(const sd_fft_ctx_t Q,
                            double* a, const double* b, ulong m_, ulong depth);
void sd_fft_ctx_point_sqr(const sd_fft_ctx_t Q,
                            double* a, ulong m_, ulong depth);


/*
    The bit reversed table is
        w = {e(0), e(1/2), e(1/4), e(3/4), e(1/8), e(5/8), e(3/8), e(7/8), ...}
    Only the terms of even index are explicitly stored, and they are split
    among several tables.
*/

/* look up w[2*j] */
FLINT_FORCE_INLINE double sd_fft_ctx_w2(const sd_fft_ctx_t Q, ulong j)
{
    ulong j_bits, j_r;
    SET_J_BITS_AND_J_R(j_bits, j_r, j);
    return Q->w2tab[j_bits][j_r];
}

/* look up -w[2*j]^-1 */
FLINT_FORCE_INLINE double sd_fft_ctx_w2inv(const sd_fft_ctx_t Q, ulong j)
{
    ulong j_bits, j_mr;
    SET_J_BITS_AND_J_MR(j_bits, j_mr, j);
    return (j == 0) ? -1.0 : Q->w2tab[j_bits][j_mr];
}

/* look up w[j] */
FLINT_FORCE_INLINE double sd_fft_ctx_w(const sd_fft_ctx_t Q, ulong j)
{
    double r = sd_fft_ctx_w2(Q, j/2);
    return (j&1) ? -r : r;
}

typedef struct {
    ulong prime;
    ulong coeff_len;
    ulong nprimes;
    ulong* data;
} crt_data_struct;

typedef crt_data_struct crt_data_t[1];

void crt_data_init(crt_data_t C, ulong prime, ulong coeff_len, ulong nprimes);

void crt_data_clear(crt_data_t C);

/* return mpn of length C->coeff_len */
FLINT_FORCE_INLINE ulong* crt_data_co_prime(const crt_data_t C, ulong i)
{
    FLINT_ASSERT(i < C->nprimes);
    return C->data + i*C->coeff_len;
}

FLINT_FORCE_INLINE ulong* _crt_data_co_prime(const crt_data_t C, ulong i, ulong n)
{
    FLINT_ASSERT(i < C->nprimes);
    FLINT_ASSERT(n == C->coeff_len);
    return C->data + i*n;
}

/* return mpn of length C->coeff_len */
FLINT_FORCE_INLINE ulong* crt_data_prod_primes(const crt_data_t C)
{
    return C->data + C->nprimes*C->coeff_len;
}

/* the reduction of co_prime mod the i^th prime */
FLINT_FORCE_INLINE ulong* crt_data_co_prime_red(const crt_data_t C, ulong i)
{
    FLINT_ASSERT(i < C->nprimes);
    return C->data + C->nprimes*C->coeff_len + C->coeff_len + i;
}


typedef void (*to_ffts_func)(
    sd_fft_ctx_struct* Qffts, double* d, ulong dstride,
    const ulong* a_, ulong an_, ulong atrunc,
    const vec4d* two_pow,
    ulong start_easy, ulong stop_easy,
    ulong start_hard, ulong stop_hard);

typedef struct {
    ulong np;
    ulong bits;
    ulong bn_bound;
    to_ffts_func to_ffts;
} profile_entry_struct;

typedef profile_entry_struct profile_entry_t[1];

#define MPN_CTX_NCRTS 8
#define MAX_NPROFILES 20
#define VEC_SZ 4

/*
    The tables for powers of two each have this fixed length. This has to go up
    linearly with the max number of primes MPN_CTX_NCRTS involved in chinese
    remaindering. This length is checked with asserts in the code.
*/
#define MPN_CTX_TWO_POWER_TAB_SIZE 256

typedef struct {
    sd_fft_ctx_struct ffts[MPN_CTX_NCRTS];
    crt_data_struct crts[MPN_CTX_NCRTS];

    /*
        For each table of tables of powers of two, the whole collection is held
        in one big buffer and the table is an array of pointer into it.
    */
    vec4d* vec_two_pow_tab[(MPN_CTX_NCRTS + VEC_SZ - 1)/VEC_SZ];
    vec4d* vec_two_pow_buffer;
    double* slow_two_pow_tab[MPN_CTX_NCRTS];
    double* slow_two_pow_buffer;

    profile_entry_struct profiles[MAX_NPROFILES];
    ulong profiles_size;
    void* buffer;
    ulong buffer_alloc;
} mpn_ctx_struct;

typedef mpn_ctx_struct mpn_ctx_t[1];

void _convert_block(ulong* Xs, sd_fft_ctx_struct* Rffts, double* d, ulong dstride, ulong np, ulong Iv);
ulong flint_mpn_nbits(const ulong* a, ulong an);
int flint_mpn_cmp_ui_2exp(const ulong* a, ulong an, ulong b, ulong e);
unsigned char flint_mpn_add_inplace_c(ulong* z, ulong zn, ulong* a, ulong an, unsigned char cf);


void mpn_ctx_init(mpn_ctx_t R, ulong p);
void mpn_ctx_clear(mpn_ctx_t R);
void* mpn_ctx_fit_buffer(mpn_ctx_t R, ulong n);
void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn);

void _nmod_poly_mul_mid_mpn_ctx(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R);

void _nmod_poly_divrem_mpn_ctx(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R);

void _nmod_poly_mul_mid_classical(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod);

void _nmod_poly_mul_mid(
    ulong* z, slong zl, slong zh,
    const ulong* a, slong an,
    const ulong* b, slong bn,
    nmod_t mod);

typedef struct {
    ulong depth;
    ulong N;
    ulong offset;
    ulong np;
    ulong stride;
    ulong bn;
    ulong btrunc;
    double* bbuf;
} mul_precomp_struct;

void _mul_precomp_init(
    mul_precomp_struct* M,
    const ulong * b, ulong bn, ulong btrunc,
    ulong depth,
    nmod_t mod,
    mpn_ctx_t R);

FLINT_FORCE_INLINE void _mul_precomp_clear(mul_precomp_struct* M)
{
    flint_aligned_free(M->bbuf);
}

int _nmod_poly_mul_mid_precomp(
    ulong* z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    mul_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R);

typedef struct {
    mul_precomp_struct quo_maker[1];
    mul_precomp_struct rem_maker[1];
} nmod_poly_divrem_precomp_struct;

FLINT_FORCE_INLINE void _nmod_poly_divrem_precomp_clear(nmod_poly_divrem_precomp_struct* M)
{
    _mul_precomp_clear(M->quo_maker);
    _mul_precomp_clear(M->rem_maker);
}

void _nmod_poly_divrem_precomp_init(
    nmod_poly_divrem_precomp_struct* M,
    const ulong* b, ulong bn,
    ulong Bn,
    nmod_t mod,
    mpn_ctx_t R);


int _nmod_poly_divrem_precomp(
    ulong* q,
    ulong* r,
    const ulong* a, ulong an,
    nmod_poly_divrem_precomp_struct* M,
    nmod_t mod,
    mpn_ctx_t R);

mpn_ctx_struct * get_default_mpn_ctx(void);

void mpn_mul_default_mpn_ctx(nn_ptr r1, nn_srcptr i1, slong n1, nn_srcptr i2, slong n2);
void _nmod_poly_mul_mid_default_mpn_ctx(nn_ptr res, slong zl, slong zh, nn_srcptr a, slong an, nn_srcptr b, slong bn, nmod_t mod);


int _fmpz_poly_mul_mid_mpn_ctx(
    fmpz * z, ulong zl, ulong zh,
    const fmpz * a, ulong an,
    const fmpz * b, ulong bn,
    mpn_ctx_t R);

int _fmpz_poly_mul_mid_default_mpn_ctx(
    fmpz * z, slong zl, slong zh,
    const fmpz * a, slong an,
    const fmpz * b, slong bn);

#ifdef __cplusplus
}
#endif

#endif

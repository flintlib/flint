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

/* alignment for transform data buffers */
#define FLINT_FFT_SMALL_ALIGNMENT 4096

/* Upper bound, in bytes, on the transform storage a single transformed
   ring or context may plan for (its expected number of simultaneously
   live elements times the per-element data size). Constructors decline
   above the bound so drivers fall back to slower algorithms rather than
   exhaust memory; mutable for tuning. A future FLINT-wide project may
   replace this with allocation-failure interception. */
FLINT_DLL extern ulong flint_fft_small_max_transformed_ring_size;

#if FLINT_HAVE_FFT_SMALL

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

FLINT_FORCE_INLINE ulong n_nbits(ulong x) {
    return FLINT_BIT_COUNT(x);
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
    Let `e(theta) = exp(2*pi*i*theta)`, let `++` denotes concatenation, and `*`
    denotes pointwise multiplication.
    Consider an infinite sequence `w` of complex numbers defined as follows.

        w_0 = {1}
        w_i = w_{i-1} ++ e(2^-i) * w_{i-1}

    For example

        w_1 = {1, -1}
        w_2 = {1, -1, i, -i}
        w_3 = {1, -1, i, -i, e(1/8), e(5/8), e(3/8), e(7/8)}

    Note that the `2i+1`-th element of `w` is the negation of the `2i`-th element,
    so we only need to store every other element of `w`. Specific elements of `w`
    can be directly accessed with `sd_fft_ctx_w`.

    In actual implementation, we work modulo a prime `p` such that `p-1` has
    high, but finite, 2-valuation. Most of the calculation works there as well,
    except the sequence `w` cannot be infinite.

    The twiddle factors are split across SD_FFT_CTX_W2TAB_SIZE tables:

        [0] = {e(0)}                                original index 0
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
    Up to the first SD_FFT_CTX_W2TAB_INIT tables are stored consecutively to ease
    the lookup of small indices, which must currently be at least max(4,
    LG_BLK_SZ).
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

/*
    This context is the one expected to sit in a global position.
    One container of this context is the thread_local storage for mpn_ctx_t.
    See sd_fft_ctx_init_prime for the meaning of each field.
*/
typedef struct {
    double p;   /* the FFT prime */
    double pinv;
    nmod_t mod;
    ulong primitive_2power_root; /* primitive 2^v-th root, v = valuation(p - 1, 2) */
#if FLINT_USES_PTHREAD
    _Atomic(unsigned int) w2tab_depth;
#else
    unsigned int w2tab_depth;
#endif
    double* w2tab[SD_FFT_CTX_W2TAB_SIZE];
    /* see above for the layout of w2tab. Each entry is an integer
       in the range [-p/2, p/2]. */
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

/* Return the address of block `I` in `d`, where each block has `BLK_SZ` entries. */
FLINT_FORCE_INLINE double* sd_fft_ctx_blk_index(double* d, ulong I)
{
    return d + sd_fft_ctx_blk_offset(I);
}

/*
location of the bit-reversed eval:
    with out_data = fft(in_data) of length 2^L, then
    eval_poly(in_data, sd_fft_ctx_w(Q, i)) = out_data[sd_fft_ctx_trunc_index(L, i)]
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
    The recursive definition of the infinite sequence `w` can be found above.
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

/* look up w[j]. Definition of w above */
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
    /* reservation of a stable tail of the scratch buffer (see
       mpn_ctx_fit_buffer_reserve); zero when none is live */
    ulong reserved_head;
    ulong reserved_tail;
    ulong buffer_alloc;

    /* constants of the two-prime pipeline: inv(p1) mod p2 with 2^-d
       folded in per depth, 2^-d mod p1, and the exact bit size of
       p1*p2 for the chunk-size bound */
    ulong np2_s1[32];       /* 2^-d mod p1 */
    ulong np2_s2[32];       /* 2^-d mod p2 */
    ulong np2_ip1;          /* inv(p1) mod p2 */
    ulong np2_prodbits;     /* floor(log2(p1*p2)) */
} mpn_ctx_struct;

typedef mpn_ctx_struct mpn_ctx_t[1];

void _convert_block(ulong* Xs, sd_fft_ctx_struct* Rffts, double* d, ulong dstride, ulong np, ulong Iv);
void _convert_block_1_mod(ulong* Xs, sd_fft_ctx_struct* Rffts, double* d, ulong dstride, ulong I, double p2, double p2inv);
ulong flint_mpn_nbits(const ulong* a, ulong an);
int flint_mpn_cmp_ui_2exp(const ulong* a, ulong an, ulong b, ulong e);
unsigned char flint_mpn_add_inplace_c(ulong* z, ulong zn, ulong* a, ulong an, unsigned char cf);


void mpn_ctx_init(mpn_ctx_t R, ulong p);
void mpn_ctx_clear(mpn_ctx_t R);
void* mpn_ctx_fit_buffer(mpn_ctx_t R, ulong n);
/* Reserve a stable tail region of the scratch buffer: the buffer is
   grown once to hold head + tail bytes and the returned tail region
   then stays at a fixed address for as long as the reservation is
   live, provided every interleaved mpn_ctx_fit_buffer request stays
   within head bytes - the caller's exclusivity contract, asserted in
   debug builds. One reservation at a time. */
void * mpn_ctx_fit_buffer_reserve(mpn_ctx_t R, ulong head, ulong tail);
void mpn_ctx_fit_buffer_release(mpn_ctx_t R);
void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong* z, const ulong* a, ulong an, const ulong* b, ulong bn);
/* Size range in which the two-prime pipeline is preferred: it beats
   the wide pipeline at every measured size from 16 limbs (where the
   wide conversions cost 2.2x while the two-prime transforms sit on
   their depth-8 padding floor) through 12000; the lower bound is 1
   since nothing below 16 changes structurally and the test suite
   covers single-limb operands, and the upper bound is measured (no
   crossover back to the wide geometry was found up to 12000, with
   the exactness bound enforcing its own independent limit). These
   are the single home for the constants used by the multiplication
   dispatch, the range dispatch and the plan gate. */
#define FLINT_MPN_MUL_NP2_MIN_BN 1
#define FLINT_MPN_MUL_NP2_MAX_AN 12000

/* two-prime pipeline for the Toom boundary; an >= bn >= 1 */
void _mpn_ctx_mpn_mul_np2(mpn_ctx_t R, nn_ptr z, nn_srcptr a, ulong an,
                    nn_srcptr b, ulong bn);
/* low zh limbs of a*b, exact */
void _mpn_ctx_mpn_mullow_np2(mpn_ctx_t R, nn_ptr z, ulong zh,
                    nn_srcptr a, ulong an, nn_srcptr b, ulong bn);
/* limbs [zl, zh): exact for zl = 0, lower approximation otherwise */
void _mpn_ctx_mpn_mul_window_np2(mpn_ctx_t R, nn_ptr z, ulong zl,
                    ulong zh, nn_srcptr a, ulong an, nn_srcptr b, ulong bn);
/* Garner + register-window recomposition of two-prime lanes; see the
   definition for the contract. depth == UWORD_MAX means the lanes are
   already plain residues. */
/* destructive on d0, d1 */
void _fft_small_np2_crt_recompose(const mpn_ctx_struct * R, nn_ptr z,
                    ulong zl, ulong zh, double * d0, double * d1,
                    ulong nchunks, ulong navail8, ulong bits,
                    ulong depth);
int _fft_small_np2_choose(const mpn_ctx_struct * R, ulong an, ulong bn,
                    ulong len_bound, ulong * bits, ulong * depth);
void _fft_small_np2_pack(double * d, double * d2, nn_srcptr x, ulong xn,
                    ulong len, ulong itrunc, ulong bits);

void _mpn_ctx_mpn_mul_range(mpn_ctx_t R, ulong* z, ulong lo, ulong hi, const ulong* a, ulong an, const ulong* b, ulong bn);

/* ------------------------------------------------------------------------ */
/* FFT plans and transformed operands                                       */
/* ------------------------------------------------------------------------ */

/*
    A plan fixes everything two operands must agree on in order to be
    combined pointwise in transform space: the set of CRT primes
    (np primes starting at index offset, or a direct single-prime FFT
    modulo a user prime), the transform depth 2^depth, and the coefficient
    bound the parameters were selected for. It also carries the output
    window [zl, zh) of a length zn convolution for the operation currently
    being performed, where ztrunc is the fft truncation length (equal to
    2^depth when the wraparound trick computes a cyclic convolution).

    The output coefficients (including any pointwise accumulation done in
    transform space) must be bounded by bound_c * 2^bound_e; np is selected
    such that the product of the primes exceeds this bound.

    Forward transforms are unnormalized; the normalization factors
    m[i] = (cop_i * 2^depth)^-1 mod p_i are folded in exactly once per
    pointwise product, as sd_fft_ctx_point_mul does today.
*/
typedef struct
{
    mpn_ctx_struct * R;

    /* CRT configuration */
    ulong np;
    ulong offset;
    sd_fft_ctx_struct * ffts;   /* R->ffts, or direct_fft below */
    crt_data_struct * crts;     /* R->crts */
    sd_fft_ctx_struct direct_fft[1];
    int use_direct_fft;

    /* transform geometry */
    ulong depth;
    ulong stride;               /* n_round_up(2^depth, 128) */

    /* window of the current operation */
    ulong zl;
    ulong zh;
    ulong zn;
    ulong ztrunc;

    /* output coefficient bound bound_c * 2^bound_e */
    ulong bound_c;
    ulong bound_e;
    ulong bits;     /* mpn plans: chunk size of the bit packing; 0 otherwise */
    int sign;

    /* normalization factors, indexed 0 <= i < np */
    ulong m[MPN_CTX_NCRTS];
} fft_small_plan_struct;

typedef fft_small_plan_struct fft_small_plan_t[1];

/* Select np and offset such that prod_primes(np) >= c * 2^e, considering
   only the primes ffts[i] for i < np_max. Returns 0 on failure. */
int _fft_small_plan_set_bound(fft_small_plan_t P, ulong c, ulong e, ulong np_max);

/* Choose depth and ztrunc for the output window [zl, zh) of a length zn
   convolution whose operand truncation lengths are at most xtrunc_max. */
/* min_depth: LG_BLK_SZ for plans whose conversions are block
   structured (the chunked mpn import and export), 4 where the
   conversions handle short transforms (the nmod paths) */
void _fft_small_plan_set_window(fft_small_plan_t P,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max, ulong min_depth);

/* Same, but with the depth of P fixed. Returns 0 if the window needs a
   transform longer than 2^depth. */
int _fft_small_plan_fit_window(fft_small_plan_t P,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max, ulong min_depth);

/* Compute stride and the normalization factors m[]. Requires the CRT
   configuration and the depth to be set. */
void _fft_small_plan_set_normalizers(fft_small_plan_t P);

/*
    Plan for computing the window [zl, zh) of a convolution of length zn
    over Z/mod.n, where at most len_bound coefficient products, each of
    absolute value less than 2^prod_bits, accumulate onto a single output
    coefficient before reduction mod mod.n. For an ordinary multiplication
    of lengths an and bn this means len_bound = min(an, bn) (any upper
    bound such as bn works) and prod_bits = 2*modbits; a bilinear
    accumulation of K products computed in transform space (e.g. matrix
    multiplication with inner dimension K) passes len_bound = K*min(an, bn),
    and unreduced inputs can be accounted for through prod_bits.

    xtrunc_max is an upper bound for the operand truncation lengths
    (n_round_up(max operand length, BLK_SZ)) and direct_len gates the
    direct FFT branch (pass the length that should be compared against the
    direct FFT cutoff, or 0 to disable the branch, e.g. when the plan will
    be reused enough that the caller wants to force it with
    fft_small_plan_force_direct_fft below).
*/
int fft_small_plan_init_nmod(fft_small_plan_t P, mpn_ctx_t R,
                    ulong zl, ulong zh, ulong zn, ulong xtrunc_max,
                    ulong len_bound, ulong prod_bits, nmod_t mod,
                    ulong direct_len);

/* Plan for cyclic convolutions modulo x^N - 1, N = 2^depth, writing the
   output coefficients [0, zh). len_bound/prod_bits as above (a full
   cyclic convolution accumulates up to N products, so an ordinary
   multiplication passes len_bound = N). */
/* Plan for transforms of nonnegative integers of at most an_max and
   bn_max limbs, with prod_of_primes large enough for len_bound
   accumulated products in transform space. Returns 0 if no admissible
   prime count and chunk size exist. */
/* len_bound counts accumulated elementary products; is_signed adds
   the spare factor of two the centered signed exports require, after
   the pipeline dispatch has keyed on the accumulation count alone,
   so a two-term signed dot product (complex multiplication, 2 x 2
   matrix multiplication) reaches the two-prime plans. */
int fft_small_plan_init_mpn(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound,
                    int is_signed);

/* plan constructor with a selectable prime regime, for profiling and
   test code: np = 0 makes the automatic choice (the same as
   fft_small_plan_init_mpn), np = 2 forces the two-prime geometry
   (failing when its exactness bound or size range cannot be met),
   and any other value forces the wide selection */
int fft_small_plan_init_mpn_np(fft_small_plan_t P, mpn_ctx_t R,
                    ulong an_max, ulong bn_max, ulong len_bound,
                    int is_signed, ulong np);

int fft_small_plan_init_nmod_cyclic(fft_small_plan_t P, mpn_ctx_t R,
                    ulong depth, ulong zh,
                    ulong len_bound, ulong prod_bits, nmod_t mod);

/* Plan for products of nonnegative integers of at most *nn limbs reduced
   mod 2^(FLINT_BITS * *nn) - 1, rounding *nn up to the nearest length
   admitting an exact chunk splitting (written back). Returns 0 if none
   exists. fft_small_export_mpn yields a nonnegative representative;
   see plan.c for its size bound and the caller-side fold. */
int fft_small_plan_init_mpn_cyclic(fft_small_plan_t P, mpn_ctx_t R,
                    ulong * nn, ulong len_bound);

void fft_small_plan_clear(fft_small_plan_t P);

/*
    A transformed operand: the forward transforms of one operand modulo
    each of the np primes of a plan, stored in the same layout as the
    scratch buffers of the fused multiplication drivers (np per-prime
    buffers of stride doubles each, in sd_fft block order).

    An operand is only compatible with plans agreeing with the
    np/offset/depth/stride fields it was created with, which is checked
    with asserts. domain distinguishes unnormalized forward transforms
    (PRIMAL) from pointwise products carrying one normalization factor
    (PRODUCT); the two must not be mixed linearly.
*/

#define FFT_SMALL_OP_PRIMAL 0
#define FFT_SMALL_OP_PRODUCT 1

typedef struct
{
    double * data;
    ulong stride;
    ulong np;
    ulong offset;
    ulong depth;
    ulong itrunc;
    int domain;
    int owns_data;
} fft_small_op_struct;

typedef fft_small_op_struct fft_small_op_t[1];

void fft_small_op_init(fft_small_op_t X, const fft_small_plan_t P);
void fft_small_op_init_borrowed(fft_small_op_t X, const fft_small_plan_t P,
                    double * data);
ulong fft_small_op_sizeof_data(const fft_small_plan_t P);
void fft_small_op_clear(fft_small_op_t X);

/* pointwise operations; Z may alias A (and B where it makes sense) */
void fft_small_op_mul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P);
void fft_small_op_sqr(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_plan_t P);
void fft_small_op_addmul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P);
void fft_small_op_submul(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P);
void fft_small_op_add(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P);
void fft_small_op_sub(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_op_t B, const fft_small_plan_t P);
void fft_small_op_neg(fft_small_op_t Z, const fft_small_op_t A,
                    const fft_small_plan_t P);

/* inverse transform in place (normalization was already applied by the
   pointwise stage) */
void fft_small_ifft(fft_small_op_t X, const fft_small_plan_t P);

/* forward transform of (a, an) reading it zero-padded to length itrunc,
   with fft truncation P->ztrunc */
void fft_small_fft_nmod(fft_small_op_t X, const ulong* a, ulong an,
                    ulong itrunc, nmod_t mod, const fft_small_plan_t P);

/* chinese remaindering of an inverse-transformed operand into the output
   window [P->zl, P->zh) */
void fft_small_fft_mpn(fft_small_op_t X, const ulong* a, ulong an,
                    const fft_small_plan_t P);
/* chinese remaindering restricted to the window [zl, zh), which must lie
   inside [P->zl, P->zh); z receives zh - zl coefficients */
void fft_small_export_nmod_range(ulong* z, const fft_small_op_t X,
                    ulong zl, ulong zh, nmod_t mod, const fft_small_plan_t P);

/* Write the low zn limbs of the accumulated integer to z. The true value
   is reduced mod 2^(FLINT_BITS*zn) if it does not fit. */
/* The signed exports reconstruct slot values centered to
   (-P/2, P/2] for the product P of the plan's primes: callers must
   keep true slot values inside that range, one spare factor of two
   against the unsigned exactness bound, e.g. by doubling len_bound
   at plan creation. */
void fft_small_export_mpn_signed(ulong* z, ulong zn, int * sign,
                    const fft_small_op_t X, ulong nslots,
                    const fft_small_plan_t P);
void fft_small_export_mpn_signed_trunc(ulong* z, ulong zn, int * sign,
                    const fft_small_op_t X, ulong nslots, ulong lo_limbs,
                    const fft_small_plan_t P);
void fft_small_export_mpn(ulong* z, ulong zn, const fft_small_op_t X,
                    const fft_small_plan_t P);

void fft_small_export_nmod(ulong* z, const fft_small_op_t X, nmod_t mod,
                    const fft_small_plan_t P);

/* ------------------------------------------------------------------------ */
/* Generic two-stage convolution engine                                     */
/* ------------------------------------------------------------------------ */

/*
    The engine runs the fused per-prime pipeline shared by the fft_small
    multiplication drivers: for each prime, convert and forward-transform
    the operands, multiply pointwise with the normalization factor folded
    in, and inverse-transform (stage 1, parallelized over the primes, with
    an optional nested worker transforming b concurrently with a); then
    chinese remainder blocks of output coefficients into the destination
    (stage 2, parallelized over output ranges).

    The coefficient-type specific parts enter through mod_fn (conversion
    of an operand into a per-prime double buffer), crt_fn (chinese
    remaindering of an output range, called once per range so that np can
    be specialized at compile time by the callers), and a tuning table of
    threshold predicates, so that each type keeps its separately tuned
    thresholds.
*/

typedef void (*fft_small_mod_func)(double* xbuf, ulong xtrunc,
    const void* x, ulong xn, slong xaux,
    const sd_fft_ctx_struct* fft, const void* params);

/* Per-range output slots for crt_fn (e.g. carries of a digit
   representation), consumed by an optional serial pass after stage 2. */
#define FFT_SMALL_CRT_LOCAL 4

typedef struct
{
    ulong stop_zi;
    ulong local[FFT_SMALL_CRT_LOCAL];
} fft_small_crt_range_struct;

typedef void (*fft_small_crt_func)(void* z, ulong zl,
    ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts, ulong min_an_bn, ulong* local,
    const void* params);

typedef void (*fft_small_s2_finish_func)(void* z, ulong zl, ulong zh,
    fft_small_crt_range_struct* ranges, ulong nranges, const void* params);

typedef struct
{
    ulong (*s1_threads)(ulong np, ulong tune_n);
    int (*s1_b_worker)(ulong np, ulong tune_n, int squaring);
    int (*s2_rethread)(ulong np, ulong zn);   /* may be NULL */
} fft_small_conv_tuning;

typedef struct
{
    const void* a; ulong an; slong aaux; ulong atrunc;
    const void* b; ulong bn; slong baux; ulong btrunc;
    const double* bfft;         /* if != NULL, b is already transformed
                                   (np buffers of stride bfft_stride) and
                                   b/bn/baux/btrunc are ignored */
    ulong bfft_stride;
    int squaring;
    ulong tune_n;               /* size fed to the tuning predicates */
    ulong min_an_bn;            /* bound on the number of products
                                   accumulating on one output coefficient */
    fft_small_mod_func mod_fn;
    fft_small_crt_func crt_fn;
    fft_small_s2_finish_func s2_finish;   /* may be NULL; called serially
                                             after stage 2 with the range
                                             boundaries and crt_fn's
                                             per-range local outputs */
    const void* params;
    const fft_small_conv_tuning* tuning;
} fft_small_conv_arg_struct;

void _fft_small_conv(void* z, const fft_small_plan_t P,
                    const fft_small_conv_arg_struct* A);

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
    fft_small_plan_struct P[1];
    fft_small_op_struct bfft[1];
    ulong bn;
    ulong btrunc;
} mul_precomp_struct;

void _mul_precomp_init(
    mul_precomp_struct* M,
    const ulong * b, ulong bn, ulong btrunc,
    ulong depth,
    nmod_t mod,
    mpn_ctx_t R);

FLINT_FORCE_INLINE void _mul_precomp_clear(mul_precomp_struct* M)
{
    fft_small_op_clear(M->bfft);
    fft_small_plan_clear(M->P);
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


/* negacyclic pointwise engine (negacyclic.c) */
typedef void (*from_ffts_func)(
    ulong* z, ulong lo, ulong hi, ulong c_lo, ulong clen,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* Rcrts,
    ulong bits,
    ulong start_easy, ulong stop_easy,
    ulong* overhang, ulong* boundbuf);

to_ffts_func _mpn_ctx_to_ffts_func(ulong np, ulong bits);
from_ffts_func _mpn_ctx_from_ffts_func(ulong np);

typedef struct sd_fft_mpn_mulmod_2expp1_scratch_struct
{
    double * buf;      /* transform buffers, 2 np data blocks */
    nn_ptr w;          /* recomposition window */
    double * digs;     /* digit doubles for b <= 50, else NULL */
} sd_fft_mpn_mulmod_2expp1_scratch_struct;

typedef sd_fft_mpn_mulmod_2expp1_scratch_struct sd_fft_mpn_mulmod_2expp1_scratch_t[1];

typedef struct sd_fft_mpn_mulmod_2expp1_ctx_struct
{
    slong N, b, m, depth, np;
    nn_ptr biaspoly;     /* precomputed bias polynomial (np == 4) */
    slong whi;           /* recomposition limb count (np == 4) */
    ulong biasres[8];    /* bias mod p_k */
    ulong p[8];
    nmod_t mod[8];
    double * tw[8];      /* w^i, twisted input scale */
    double * itw[8];     /* w^-i * m^-1, output scale */
    slong * crt2;        /* scratch: two-prime corrections */
    fmpz_t Pbig, Phalf;
    mpn_ctx_struct * R;
    slong dsz;
} sd_fft_mpn_mulmod_2expp1_ctx_struct;

typedef sd_fft_mpn_mulmod_2expp1_ctx_struct sd_fft_mpn_mulmod_2expp1_ctx_t[1];

slong sd_fft_mpn_mulmod_2expp1_choose_np(mpn_ctx_struct * R, slong b, slong depth);
void sd_fft_mpn_mulmod_2expp1_ctx_init(sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            mpn_ctx_struct * R, slong N, slong m);
void sd_fft_mpn_mulmod_2expp1_ctx_clear(sd_fft_mpn_mulmod_2expp1_ctx_struct * C);
void sd_fft_mpn_mulmod_2expp1(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, nn_ptr z, nn_srcptr x,
            nn_srcptr y, sd_fft_mpn_mulmod_2expp1_scratch_struct * S);
/* staged variant: F (transformed_size doubles) caches one operand's
   transform for reuse across mul_cached calls; the cached operand
   corresponds to y of the one-shot function */
slong sd_fft_mpn_mulmod_2expp1_transformed_size(
            const sd_fft_mpn_mulmod_2expp1_ctx_struct * C);
void sd_fft_mpn_mulmod_2expp1_transform(sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            double * F, nn_srcptr y,
            sd_fft_mpn_mulmod_2expp1_scratch_struct * S);
void sd_fft_mpn_mulmod_2expp1_mul_cached(sd_fft_mpn_mulmod_2expp1_ctx_struct * C,
            nn_ptr z, nn_srcptr x, const double * F,
            sd_fft_mpn_mulmod_2expp1_scratch_struct * S);

void sd_fft_mpn_mulmod_2expp1_scratch_init(sd_fft_mpn_mulmod_2expp1_scratch_struct * S,
            const sd_fft_mpn_mulmod_2expp1_ctx_struct * C);
void sd_fft_mpn_mulmod_2expp1_scratch_clear(sd_fft_mpn_mulmod_2expp1_scratch_struct * S);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_HAVE_FFT_SMALL */

#endif

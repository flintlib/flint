/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024, 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_VEC_H
#define NMOD_VEC_H

#ifdef NMOD_VEC_INLINES_C
#define NMOD_VEC_INLINE
#else
#define NMOD_VEC_INLINE static inline
#endif

#include "flint.h"
#include "nmod.h"  /* nmod_mul, nmod_fmma */
#include "flint-mparam.h"  /* FLINT_NMOD_VEC_DOT_{U52,SPLIT_LIMBS}_MIN_LEN, ..._MAX_BITS */

/*
    SIMD dot products for moduli above 2^32 (see dot_u52.c, dot_u64.c,
    dot_split_limbs.c):
    - on AVX512-IFMA, the strategy of nmod_mat_mul_u52 for moduli up to 2^52
      (_DOT_U52, only when the unreduced dot product fits in two limbs), and
      for larger moduli the products of the 32-bit halves of the entries
      (_DOT_U64 in the two-limb band, _DOT3_U64 in the three-limb one);
    - on AVX2 / AVX-512, the entries split into two limbs of at most 32 bits
      multiplied with vpmuludq (_DOT_SPLIT_LIMBS / _DOT3_SPLIT_LIMBS), for
      moduli above 2^32 up to 2^NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS (with IFMA,
      from 53 bits on).
    The limb-count convention on dot_method_t is kept: each method appears
    in the band(s) it can serve.
*/
#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512DQ__) \
        && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_VEC_HAVE_DOT_U52 1
# define NMOD_VEC_HAVE_DOT_U64 1
#else
# define NMOD_VEC_HAVE_DOT_U52 0
# define NMOD_VEC_HAVE_DOT_U64 0
#endif

#if FLINT_BITS == 64 && defined(__AVX2__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_VEC_HAVE_DOT_SPLIT_LIMBS 1
#else
# define NMOD_VEC_HAVE_DOT_SPLIT_LIMBS 0
#endif

/*
    Minimal lengths for these to beat the scalar code, measured on Ice Lake,
    Zen 4, Cascade Lake, Broadwell and Emerald Rapids. For u64, measured
    against _DOT2 in the two-limb band and against _DOT3_ACC / _DOT3 in the
    three-limb one, the crossovers are 40-80 and 32-56 on the IFMA machines
    measured, around these compromises. The split-limbs crossover varies from 48 (Zen 4) to 128 (Cascade
    Lake) and never (Broadwell), and the u52 one from 32 (Ice Lake, Emerald
    Rapids: 1.1x and 1.5x faster than _DOT2 at lengths 32 to 39) to 40 (Zen
    4: 1.1x slower at 32), hence parameters of flint-mparam.h (for split
    limbs, 0 means never).
*/
#if !defined(FLINT_NMOD_VEC_DOT_U52_MIN_LEN) \
    || !defined(FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN) \
    || !defined(FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS)
# error "flint-mparam.h must define FLINT_NMOD_VEC_DOT_U52_MIN_LEN, FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN and FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS"
#endif
#define NMOD_VEC_DOT_U52_MIN_LEN FLINT_NMOD_VEC_DOT_U52_MIN_LEN
#define NMOD_VEC_DOT_U64_MIN_LEN 80
#define NMOD_VEC_DOT3_U64_MIN_LEN 48
#define NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN
#define NMOD_VEC_DOT_SPLIT_LIMBS_ENABLED \
    (NMOD_VEC_HAVE_DOT_SPLIT_LIMBS && NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN > 0)

/* with AVX512-IFMA, typical largest modulus size for split limbs: they
   are faster than u64 up to 60 bits on Ice Lake (61 on Zen 4, a tie at 61
   bits there); used by the tests */
#define NMOD_VEC_DOT_SPLIT_LIMBS_IFMA_MAX_BITS 60

/*
    Largest modulus bit size for split limbs, a parameter of flint-mparam.h
    (FLINT_BITS or more for no limit). Without IFMA the alternative is the
    scalar code: on Meteor Lake and Raptor Lake the split-limbs variant used
    from 62 bits on is 1.1-1.2x slower than _DOT3_ACC at every length, and
    on Meteor Lake the one of 61 bits too. With IFMA it is u64, which takes over between 59
    (Emerald Rapids) and 62 bits (Tiger Lake).
*/
#define NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS FLINT_NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS
/* largest modulus for split limbs */
#if NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS >= FLINT_BITS
# define NMOD_VEC_DOT_SPLIT_LIMBS_MAX_N UWORD_MAX
#else
# define NMOD_VEC_DOT_SPLIT_LIMBS_MAX_N (UWORD(1) << NMOD_VEC_DOT_SPLIT_LIMBS_MAX_BITS)
#endif

/*
    AArch64: the scalar _DOT2_SPLIT is faster than the two-chain _DOT2 up
    to about this length only, and _DOT2_HALF (scalar without AVX2) never
    (Apple M4); on x86 _DOT2_SPLIT and _DOT2_HALF are vectorized with AVX2
*/
#if FLINT_BITS == 64 && defined(__aarch64__)
# define NMOD_VEC_DOT2_SPLIT_MAX_LEN 160
#endif

/* below their minimal length (parameters computed for a longer length, as
   for the short dot products of a classical polynomial product), the SIMD
   methods go directly to the scalar code of their band: _DOT2 (two) or
   _DOT3_ACC / _DOT3 (three limbs) */
#if FLINT_BITS == 64
# define _NMOD_VEC_DOT3_ACC_BOUND UWORD(6521908912666391107)
#else
# define _NMOD_VEC_DOT3_ACC_BOUND UWORD(1518500250)
#endif
#define _NMOD_VEC_DOT_SCALAR(two, mod, SUF, ...)                         \
    ((two) ? _nmod_vec_dot2 ## SUF(__VA_ARGS__)                           \
     : ((mod).n <= _NMOD_VEC_DOT3_ACC_BOUND)                              \
       ? _nmod_vec_dot3_acc ## SUF(__VA_ARGS__)                           \
       : _nmod_vec_dot3 ## SUF(__VA_ARGS__))

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_VEC_NORM(vec, i)                   \
do {                                            \
    while ((i) && vec[(i) - 1] == UWORD(0))     \
        (i)--;                                  \
} while (0)

NMOD_VEC_INLINE
nn_ptr _nmod_vec_init(slong len)
{
    return (nn_ptr) flint_malloc(len * sizeof(ulong));
}

NMOD_VEC_INLINE
void _nmod_vec_clear(nn_ptr vec)
{
   flint_free(vec);
}

void _nmod_vec_randtest(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod);
void _nmod_vec_rand(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod);

NMOD_VEC_INLINE
void _nmod_vec_zero(nn_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        vec[i] = 0;
}

flint_bitcnt_t _nmod_vec_max_bits(nn_srcptr vec, slong len);

NMOD_VEC_INLINE
void _nmod_vec_set(nn_ptr res, nn_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        res[i] = vec[i];
}

NMOD_VEC_INLINE
void _nmod_vec_swap(nn_ptr a, nn_ptr b, slong length)
{
    slong i;
    for (i = 0; i < length; i++)
    {
        ulong t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}

NMOD_VEC_INLINE
int _nmod_vec_equal(nn_srcptr vec, nn_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (vec[i] != vec2[i]) return 0;

    return 1;
}

NMOD_VEC_INLINE
int _nmod_vec_is_zero(nn_srcptr vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (vec[i] != 0) return 0;

    return 1;
}

/* some IO functions */
#ifdef FLINT_HAVE_FILE
int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod);
#endif

void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod);

/* reduce, add, scalar mul */
void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod);

void _nmod_vec_add(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod);

void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_mul_nmod_redc(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_mul_nmod_generic(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_addmul_nmod_generic(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_addmul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);

#if FLINT_BITS == 64 && defined(__AVX2__)
void _nmod_vec_nored_scalar_addmul_halflimb_avx2(nn_ptr res, nn_srcptr vec, slong len, ulong c);
#elif FLINT_BITS == 64 && defined(__ARM_NEON)
void _nmod_vec_nored_scalar_addmul_halflimb_neon(nn_ptr res, nn_srcptr vec, slong len, ulong c);
#endif

NMOD_VEC_INLINE void
_nmod_vec_nored_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c)
{
#if FLINT_BITS == 64 && defined(__AVX2__)
    if (len >= 16)
    {
        _nmod_vec_nored_scalar_addmul_halflimb_avx2(res, vec, len, c);
        return;
    }
#elif FLINT_BITS == 64 && defined(__ARM_NEON)
    if (len >= 16)
    {
        _nmod_vec_nored_scalar_addmul_halflimb_neon(res, vec, len, c);
        return;
    }
#endif

    slong i;
    for (i = 0; i < len; i++)
        res[i] += vec[i] * c;
}

void _nmod_vec_nored_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c);
void _nmod_vec_nored_ll_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c);
void _nmod_vec_nored_ll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c);
void _nmod_vec_nored_lll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c);

/* invert each coefficient */
void _nmod_vec_invert(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);
void _nmod_vec_invert_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);
void _nmod_vec_invert_generic(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);
void _nmod_vec_invert_shoup(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod);

/* ---- compute dot parameters ---- */

// for _DOT2_SPLIT
#if (FLINT_BITS == 64)
#   define DOT_SPLIT_BITS 56
#   define DOT_SPLIT_MASK UWORD(72057594037927935) // (1L << DOT_SPLIT_BITS) - 1
#endif // FLINT_BITS == 64

/* The method is chosen by the library (src/nmod_vec/dot_params.c), so that
   it matches the kernels FLINT was compiled with, whatever the compiler flags
   of the calling code. */
dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod);

int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);
int _nmod_vec_dot_bound_limbs_from_params(slong len, nmod_t mod, dot_params_t params);


/* ------ dot product, specific algorithms ------ */

/* vec1[i] * vec2[i] */
ulong _nmod_vec_dot_pow2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot1(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2_half(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
ulong _nmod_vec_dot2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3_acc(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
#if FLINT_BITS == 64
ulong _nmod_vec_dot2_split(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif  // FLINT_BITS == 64
ulong _nmod_vec_dot_u52(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot_u64(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot_split_limbs(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);

/* vec1[i] * vec2[len-1-i] */
ulong _nmod_vec_dot_pow2_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot1_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2_half_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3_acc_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
#if FLINT_BITS == 64
ulong _nmod_vec_dot2_split_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif  // FLINT_BITS == 64
ulong _nmod_vec_dot_u52_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot_u64_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot_split_limbs_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);

/* vec1[i] * vec2[i][offset] */
ulong _nmod_vec_dot_pow2_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot1_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot2_half_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot2_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot3_acc_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot3_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
#if FLINT_BITS == 64
ulong _nmod_vec_dot2_split_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, ulong pow2_precomp);
#endif  // FLINT_BITS == 64
ulong _nmod_vec_dot_u52_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot_u64_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);
ulong _nmod_vec_dot_split_limbs_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod);





/* ------------- dot product, general --------------- */

// auxiliary for short dot products
// (fixedlen small constant: compiler unrolls the loops completely)
#define _NMOD_VEC_DOT_SHORT1(fixedlen,expr1,expr2)       \
        {                                                \
            ulong res = (expr1) * (expr2); i++;          \
            for (slong j = 0; j < fixedlen-1; j++, i++)  \
                res += (expr1) * (expr2);                \
            NMOD_RED(res, res, mod);                     \
            return res;                                  \
        }                                                \

#define _NMOD_VEC_DOT_SHORT2(fixedlen,expr1,expr2)          \
        {                                                   \
            ulong s0, s1, u0, u1;                           \
            umul_ppmm(u1, u0, (expr1), (expr2)); i++;       \
            for (slong j = 0; j < fixedlen-1; j++, i++)     \
            {                                               \
                umul_ppmm(s1, s0, (expr1), (expr2));        \
                add_ssaaaa(u1, u0, u1, u0, s1, s0);         \
            }                                               \
            NMOD2_RED2(s0, u1, u0, mod);                    \
            return s0;                                      \
        }                                                   \

#define _NMOD_VEC_DOT_SHORT3(fixedlen,expr1,expr2)           \
        {                                                    \
            ulong t2 = UWORD(0);                             \
            ulong t1, t0;                                    \
            umul_ppmm(t1, t0, (expr1), (expr2)); i++;        \
            for (slong j = 0; j < fixedlen - 1; j++, i++)    \
            {                                                \
                ulong s0, s1;                                \
                umul_ppmm(s1, s0, (expr1), (expr2));         \
                add_sssaaaaaa(t2, t1, t0,                    \
                                t2, t1, t0,                  \
                                UWORD(0), s1, s0);           \
            }                                                \
                                                             \
            NMOD_RED(t2, t2, mod);                           \
            ulong res;                                       \
            NMOD_RED3(res, t2, t1, t0, mod);                 \
            return res;                                      \
        }                                                    \

/* - supports 1 <= len <= 11, requires method==DOT1|DOT2|DOT3|DOT_POW2          */
/* - the SIMD methods are left out on purpose (a short length with them only    */
/*   occurs with parameters computed for a longer length, and their functions   */
/*   fall back to the scalar code): testing them here too made can slow down    */
/*   the callers' loops (e.g. _nmod_poly_mul_classical) at short lengths        */
/* - i must be already initialized at the first wanted value                    */
#define _NMOD_VEC_DOT_SHORT(i, expr1, expr2, len, mod, method)          \
{                                                                       \
    if (method == _DOT1 || method == _DOT_POW2)                         \
    {                                                                   \
        if (len ==  1) _NMOD_VEC_DOT_SHORT1( 1, expr1, expr2)           \
        if (len ==  2) _NMOD_VEC_DOT_SHORT1( 2, expr1, expr2)           \
        if (len ==  3) _NMOD_VEC_DOT_SHORT1( 3, expr1, expr2)           \
        if (len ==  4) _NMOD_VEC_DOT_SHORT1( 4, expr1, expr2)           \
        if (len ==  5) _NMOD_VEC_DOT_SHORT1( 5, expr1, expr2)           \
        if (len ==  6) _NMOD_VEC_DOT_SHORT1( 6, expr1, expr2)           \
        if (len ==  7) _NMOD_VEC_DOT_SHORT1( 7, expr1, expr2)           \
        if (len ==  8) _NMOD_VEC_DOT_SHORT1( 8, expr1, expr2)           \
        if (len ==  9) _NMOD_VEC_DOT_SHORT1( 9, expr1, expr2)           \
        if (len == 10) _NMOD_VEC_DOT_SHORT1(10, expr1, expr2)           \
        _NMOD_VEC_DOT_SHORT1(11, expr1, expr2)                          \
    }                                                                   \
                                                                        \
    else if (method == _DOT2)                                           \
    {                                                                   \
        if (len ==  1) return nmod_mul((expr1), (expr2), mod);          \
        if (len ==  2) _NMOD_VEC_DOT_SHORT2( 2, expr1, expr2)           \
        if (len ==  3) _NMOD_VEC_DOT_SHORT2( 3, expr1, expr2)           \
        if (len ==  4) _NMOD_VEC_DOT_SHORT2( 4, expr1, expr2)           \
        if (len ==  5) _NMOD_VEC_DOT_SHORT2( 5, expr1, expr2)           \
        if (len ==  6) _NMOD_VEC_DOT_SHORT2( 6, expr1, expr2)           \
        if (len ==  7) _NMOD_VEC_DOT_SHORT2( 7, expr1, expr2)           \
        if (len ==  8) _NMOD_VEC_DOT_SHORT2( 8, expr1, expr2)           \
        if (len ==  9) _NMOD_VEC_DOT_SHORT2( 9, expr1, expr2)           \
        if (len == 10) _NMOD_VEC_DOT_SHORT2(10, expr1, expr2)           \
        _NMOD_VEC_DOT_SHORT2(11, expr1, expr2)                          \
    }                                                                   \
                                                                        \
    else if (method == _DOT3)                                           \
    {                                                                   \
        if (len ==  1) return nmod_mul((expr1), (expr2), mod);          \
        if (len ==  2) _NMOD_VEC_DOT_SHORT3( 2, expr1, expr2)           \
        if (len ==  3) _NMOD_VEC_DOT_SHORT3( 3, expr1, expr2)           \
        if (len ==  4) _NMOD_VEC_DOT_SHORT3( 4, expr1, expr2)           \
        if (len ==  5) _NMOD_VEC_DOT_SHORT3( 5, expr1, expr2)           \
        if (len ==  6) _NMOD_VEC_DOT_SHORT3( 6, expr1, expr2)           \
        if (len ==  7) _NMOD_VEC_DOT_SHORT3( 7, expr1, expr2)           \
        if (len ==  8) _NMOD_VEC_DOT_SHORT3( 8, expr1, expr2)           \
        if (len ==  9) _NMOD_VEC_DOT_SHORT3( 9, expr1, expr2)           \
        if (len == 10) _NMOD_VEC_DOT_SHORT3(10, expr1, expr2)           \
        _NMOD_VEC_DOT_SHORT3(11, expr1, expr2)                          \
    }                                                                   \
}  while(0);                                                            \

FLINT_FORCE_INLINE ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    if (len <= 11)
    {
        if (len == 0) return UWORD(0);
        slong i = 0;
        _NMOD_VEC_DOT_SHORT(i, vec1[i], vec2[i], len, mod, params.method);
    }

    if (params.method == _DOT1)
        return _nmod_vec_dot1(vec1, vec2, len, mod);

#if FLINT_BITS == 64
    if (params.method == _DOT2_SPLIT)
        return _nmod_vec_dot2_split(vec1, vec2, len, mod, params.pow2_precomp);
#endif // FLINT_BITS == 64

    if (params.method == _DOT2)
        return _nmod_vec_dot2(vec1, vec2, len, mod);

    if (params.method == _DOT_U52)
        return (len >= NMOD_VEC_DOT_U52_MIN_LEN)
            ? _nmod_vec_dot_u52(vec1, vec2, len, mod)
            : _nmod_vec_dot2(vec1, vec2, len, mod);

    if (params.method == _DOT_SPLIT_LIMBS || params.method == _DOT3_SPLIT_LIMBS)
        return (len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
            ? _nmod_vec_dot_split_limbs(vec1, vec2, len, mod)
            : _NMOD_VEC_DOT_SCALAR(params.method == _DOT_SPLIT_LIMBS, mod, , vec1, vec2, len, mod);

    if (params.method == _DOT_U64 || params.method == _DOT3_U64)
        return (len >= ((params.method == _DOT_U64) ? NMOD_VEC_DOT_U64_MIN_LEN
                                                    : NMOD_VEC_DOT3_U64_MIN_LEN))
            ? _nmod_vec_dot_u64(vec1, vec2, len, mod)
            : _NMOD_VEC_DOT_SCALAR(params.method == _DOT_U64, mod, , vec1, vec2, len, mod);

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc(vec1, vec2, len, mod);

    if (params.method == _DOT3)
        return _nmod_vec_dot3(vec1, vec2, len, mod);

    if (params.method == _DOT2_HALF)
        return _nmod_vec_dot2_half(vec1, vec2, len, mod, params.pow2_precomp);

    if (params.method == _DOT_POW2)
    {
        if (mod.n <= UWORD(1) << (FLINT_BITS / 2))
            return _nmod_vec_dot1(vec1, vec2, len, mod);
        else
            return _nmod_vec_dot_pow2(vec1, vec2, len, mod);
    }

    // covers _DOT0 for len > 11 (i.e. mod.n == 1...)
    return UWORD(0);
}

FLINT_FORCE_INLINE ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    if (len <= 11)
    {
        if (len == 0) return UWORD(0);
        slong i = 0;
        _NMOD_VEC_DOT_SHORT(i, vec1[i], vec2[len-1-i], len, mod, params.method);
    }

    if (params.method == _DOT1)
        return _nmod_vec_dot1_rev(vec1, vec2, len, mod);

#if FLINT_BITS == 64
    if (params.method == _DOT2_SPLIT)
        return _nmod_vec_dot2_split_rev(vec1, vec2, len, mod, params.pow2_precomp);
#endif // FLINT_BITS == 64

    if (params.method == _DOT2)
        return _nmod_vec_dot2_rev(vec1, vec2, len, mod);

    if (params.method == _DOT_U52)
        return (len >= NMOD_VEC_DOT_U52_MIN_LEN)
            ? _nmod_vec_dot_u52_rev(vec1, vec2, len, mod)
            : _nmod_vec_dot2_rev(vec1, vec2, len, mod);

    if (params.method == _DOT_SPLIT_LIMBS || params.method == _DOT3_SPLIT_LIMBS)
        return (len >= NMOD_VEC_DOT_SPLIT_LIMBS_MIN_LEN)
            ? _nmod_vec_dot_split_limbs_rev(vec1, vec2, len, mod)
            : _NMOD_VEC_DOT_SCALAR(params.method == _DOT_SPLIT_LIMBS, mod, _rev, vec1, vec2, len, mod);

    if (params.method == _DOT_U64 || params.method == _DOT3_U64)
        return (len >= ((params.method == _DOT_U64) ? NMOD_VEC_DOT_U64_MIN_LEN
                                                    : NMOD_VEC_DOT3_U64_MIN_LEN))
            ? _nmod_vec_dot_u64_rev(vec1, vec2, len, mod)
            : _NMOD_VEC_DOT_SCALAR(params.method == _DOT_U64, mod, _rev, vec1, vec2, len, mod);

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc_rev(vec1, vec2, len, mod);

    if (params.method == _DOT3)
        return _nmod_vec_dot3_rev(vec1, vec2, len, mod);

    if (params.method == _DOT2_HALF)
        return _nmod_vec_dot2_half_rev(vec1, vec2, len, mod);

    if (params.method == _DOT_POW2)
    {
        if (mod.n <= UWORD(1) << (FLINT_BITS / 2))
            return _nmod_vec_dot1_rev(vec1, vec2, len, mod);
        else
            return _nmod_vec_dot_pow2_rev(vec1, vec2, len, mod);
    }

    // covers _DOT0 for len > 11 (i.e. mod.n == 1...)
    return UWORD(0);
}

FLINT_FORCE_INLINE ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, dot_params_t params)
{
    if (len <= 11)
    {
        if (len == 0) return UWORD(0);
        slong i = 0;
        _NMOD_VEC_DOT_SHORT(i, vec1[i], vec2[i][offset], len, mod, params.method);
    }

    if (params.method == _DOT1)
        return _nmod_vec_dot1_ptr(vec1, vec2, offset, len, mod);

#if FLINT_BITS == 64
    if (params.method == _DOT2_SPLIT)
        return _nmod_vec_dot2_split_ptr(vec1, vec2, offset, len, mod, params.pow2_precomp);
#endif // FLINT_BITS == 64

    /* the SIMD methods are not used here: their _ptr variants load vec2
       through gathers, slower than this scalar code on most machines
       measured (Skylake to Ice Lake, where the microcode mitigation
       of Gather Data Sampling makes them 2-3x slower; but also Zen4
       although it was not affected by this mitigation) */
    if (params.method == _DOT2 || params.method == _DOT_U52
            || params.method == _DOT_SPLIT_LIMBS || params.method == _DOT_U64)
        return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc_ptr(vec1, vec2, offset, len, mod);

    if (params.method == _DOT3 || params.method == _DOT3_U64
            || params.method == _DOT3_SPLIT_LIMBS)
    {
        if (mod.n <= _NMOD_VEC_DOT3_ACC_BOUND)  /* as in _nmod_vec_dot_params */
            return _nmod_vec_dot3_acc_ptr(vec1, vec2, offset, len, mod);
        return _nmod_vec_dot3_ptr(vec1, vec2, offset, len, mod);
    }

    if (params.method == _DOT2_HALF)
        return _nmod_vec_dot2_half_ptr(vec1, vec2, offset, len, mod);

    if (params.method == _DOT_POW2)
    {
        if (mod.n <= UWORD(1) << (FLINT_BITS / 2))
            return _nmod_vec_dot1_ptr(vec1, vec2, offset, len, mod);
        else
            return _nmod_vec_dot_pow2_ptr(vec1, vec2, offset, len, mod);
    }

    // covers _DOT0 for len > 11 (i.e. mod.n == 1...)
    return UWORD(0);
}

#undef _NMOD_VEC_DOT_SHORT1
#undef _NMOD_VEC_DOT_SHORT2
#undef _NMOD_VEC_DOT_SHORT3


/* ---- macros for dot product with expressions, specific algorithms ---- */

// _DOT1   (1 limb)
#define _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    res = UWORD(0);                                                   \
    for (i = 0; i < (len); i++)                                       \
        res += (expr1) * (expr2);                                     \
    NMOD_RED(res, res, mod);                                          \
} while(0);

// _DOT2_SPLIT   (2 limbs, splitting at DOT_SPLIT_BITS bits, 8-unrolling)
#if (FLINT_BITS == 64)
#define _NMOD_VEC_DOT2_SPLIT(res, i, len, expr1, expr2, mod, pow2_precomp) \
do                                                    \
{                                                     \
    ulong dp_lo = 0;                                  \
    ulong dp_hi = 0;                                  \
                                                      \
    for (i = 0; i+7 < (len); )                        \
    {                                                 \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
        dp_lo += (expr1) * (expr2); i++;              \
                                                      \
        dp_hi += dp_lo >> DOT_SPLIT_BITS;             \
        dp_lo &= DOT_SPLIT_MASK;                      \
    }                                                 \
                                                      \
    for ( ; i < (len); i++)                           \
        dp_lo += (expr1) * (expr2);                   \
                                                      \
    res = pow2_precomp * dp_hi + dp_lo;               \
    NMOD_RED(res, res, mod);                          \
} while(0);
#endif  // FLINT_BITS == 64

// _DOT2_HALF   (two limbs, modulus < 2**32)
// mod.n is too close to 2**32 to accumulate in some ulong
// still interesting: a bit faster than _NMOD_VEC_DOT2
#define _NMOD_VEC_DOT2_HALF(res, i, len, expr1, expr2, mod)           \
do                                                                    \
{                                                                     \
    ulong s0zz = UWORD(0);                                            \
    ulong s1zz = UWORD(0);                                            \
    for (i = 0; i < (len); i++)                                       \
    {                                                                 \
        const ulong prodzz = (expr1) * (expr2);                       \
        add_ssaaaa(s1zz, s0zz, s1zz, s0zz, 0, prodzz);                \
    }                                                                 \
    NMOD2_RED2(res, s1zz, s0zz, mod);                                 \
} while(0);

// _DOT2   (two limbs, general)
// two independent accumulators (even / odd terms): the additions with carry
// of consecutive terms do not wait for each other
#define _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    ulong u0zz = UWORD(0);                                            \
    ulong u1zz = UWORD(0);                                            \
    ulong v0zz = UWORD(0);                                            \
    ulong v1zz = UWORD(0);                                            \
                                                                      \
    for (i = 0; i+7 < (len); )                                        \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(v1zz, v0zz, v1zz, v0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(v1zz, v0zz, v1zz, v0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(v1zz, v0zz, v1zz, v0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(v1zz, v0zz, v1zz, v0zz, s1zz, s0zz);               \
        i++;                                                          \
    }                                                                 \
    for ( ; i < (len); i++)                                           \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
    }                                                                 \
                                                                      \
    add_ssaaaa(u1zz, u0zz, u1zz, u0zz, v1zz, v0zz);                   \
    NMOD2_RED2(res, u1zz, u0zz, mod);                                 \
} while(0);

// _DOT3_ACC   (three limbs, delayed accumulations)
// 8-unroll: requires  8 * (mod.n - 1)**2 < 2**128
#define _NMOD_VEC_DOT3_ACC(res, i, len, expr1, expr2, mod)            \
do                                                                    \
{                                                                     \
    ulong t2zz = UWORD(0);                                            \
    ulong t1zz = UWORD(0);                                            \
    ulong t0zz = UWORD(0);                                            \
                                                                      \
    for (i = 0; i+7 < (len); )                                        \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        ulong u0zz = UWORD(0);                                        \
        ulong u1zz = UWORD(0);                                        \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
        i++;                                                          \
        add_sssaaaaaa(t2zz, t1zz, t0zz,                               \
                        t2zz, t1zz, t0zz,                             \
                        UWORD(0), u1zz, u0zz);                        \
    }                                                                 \
                                                                      \
    ulong s0zz, s1zz;                                                 \
    ulong u0zz = UWORD(0);                                            \
    ulong u1zz = UWORD(0);                                            \
    for ( ; i < (len); i++)                                           \
    {                                                                 \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
    }                                                                 \
                                                                      \
    add_sssaaaaaa(t2zz, t1zz, t0zz,                                   \
                    t2zz, t1zz, t0zz,                                 \
                    UWORD(0), u1zz, u0zz);                            \
                                                                      \
    NMOD_RED(t2zz, t2zz, mod);                                        \
    NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                            \
} while(0);

// _DOT3   (three limbs, general)
// mod.n is too close to 2**64 to accumulate in two words
#define _NMOD_VEC_DOT3(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    ulong t2zz = UWORD(0);                                            \
    ulong t1zz = UWORD(0);                                            \
    ulong t0zz = UWORD(0);                                            \
    for (i = 0; i < (len); i++)                                       \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_sssaaaaaa(t2zz, t1zz, t0zz,                               \
                        t2zz, t1zz, t0zz,                             \
                        UWORD(0), s1zz, s0zz);                        \
    }                                                                 \
                                                                      \
    NMOD_RED(t2zz, t2zz, mod);                                        \
    NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                            \
} while(0);

/* ---- macros for dot product with expressions, general ---- */
// currently no vectorization here

#if (FLINT_BITS == 64)

#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, params)          \
do                                                                    \
{                                                                     \
    res = UWORD(0);   /* covers _DOT0 */                              \
    if (params.method == _DOT1 || params.method == _DOT_POW2)         \
        _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                \
    else if (params.method == _DOT2_SPLIT)                            \
        _NMOD_VEC_DOT2_SPLIT(res, i, len, expr1, expr2, mod,          \
                params.pow2_precomp)                                  \
    else if (params.method == _DOT2_HALF)                             \
        _NMOD_VEC_DOT2_HALF(res, i, len, expr1, expr2, mod)           \
    else if (params.method == _DOT2 || params.method == _DOT_U52      \
             || params.method == _DOT_SPLIT_LIMBS                     \
             || params.method == _DOT_U64)                            \
        _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
    else if (params.method == _DOT3_ACC                               \
             || ((params.method == _DOT3_U64                          \
                  || params.method == _DOT3_SPLIT_LIMBS)              \
                 && (mod).n <= _NMOD_VEC_DOT3_ACC_BOUND))             \
        _NMOD_VEC_DOT3_ACC(res, i, len, expr1, expr2, mod)            \
    else if (params.method == _DOT3 || params.method == _DOT3_U64     \
             || params.method == _DOT3_SPLIT_LIMBS)                   \
        _NMOD_VEC_DOT3(res, i, len, expr1, expr2, mod)                \
} while(0);

#else  // FLINT_BITS == 64

#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, params)          \
do                                                                    \
{                                                                     \
    res = UWORD(0);   /* covers _DOT0 */                              \
    if (params.method == _DOT1 || params.method == _DOT_POW2)         \
        _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                \
    else if (params.method == _DOT2_HALF)                             \
        _NMOD_VEC_DOT2_HALF(res, i, len, expr1, expr2, mod)           \
    else if (params.method == _DOT2)                                  \
        _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
    else if (params.method == _DOT3_ACC)                              \
        _NMOD_VEC_DOT3_ACC(res, i, len, expr1, expr2, mod)            \
    else if (params.method == _DOT3)                                  \
        _NMOD_VEC_DOT3(res, i, len, expr1, expr2, mod)                \
} while(0);

#endif  // FLINT_BITS == 64

/* Note: if de-inlining this, consider specializing for the common use
   case of stride1 == 1 or stride2 == 1. */
NMOD_VEC_INLINE ulong _nmod_vec_dot_strided(nn_srcptr vec1, slong stride1,
    nn_srcptr vec2, slong stride2, slong len, nmod_t mod, dot_params_t params)
{
    slong i;
    ulong r;
    NMOD_VEC_DOT(r, i, len, vec1[i * stride1], vec2[i * stride2], mod, params)
    return r;
}

#ifdef __cplusplus
}
#endif

#endif

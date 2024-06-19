/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Vincent Neiger

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

void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod);

void _nmod_vec_add(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod);

void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);
void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);

/* -------------------- dot product ----------------------- */
/* more comments in nmod_vec/dot.c */

// for _DOT2_SPLIT (currently restricted to FLINT_BITS == 64)
#if (FLINT_BITS == 64)
#   define DOT_SPLIT_BITS 56
#   define DOT_SPLIT_MASK UWORD(72057594037927935) // (1L << DOT_SPLIT_BITS) - 1
#endif // FLINT_BITS == 64

typedef enum
{
    _DOT0 = 0,           /* len == 0 || mod.n == 1 */
    _DOT_POW2 = 1,       /* modulus is a power of 2, computations performed on 1 limb */
    _DOT1 = 2,           /* 1 limb */
#if (FLINT_BITS == 64)
    _DOT2_SPLIT = 3,     /* 2 limbs, modulus < ~2**30.5 (FLINT_BITS == 64 only) */
#endif  // FLINT_BITS == 64
    _DOT2_HALF = 4,      /* 2 limbs, modulus < 2**(FLINT_BITS/2) */
    _DOT2 = 5,           /* 2 limbs */
    _DOT3_ACC = 6,       /* 3 limbs, modulus < 2**62.5 allowing accumulation in 2 limbs */
    _DOT3 = 7,           /* 3 limbs, modulus >= 2**62.5 */
} dot_method_t;
// if mod.n is not a power of 2, then
//      > 1 limb <-> method > _DOT1   |   > 2 limbs <-> method > _DOT2

typedef struct
{
    dot_method_t method;
    ulong pow2_precomp;  /* for splitting: (1L << 56) % mod.n */
} dot_params_t;

// compute dot parameters
dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod);

// _DOT1   (1 limb)
#define _NMOD_VEC_DOT1(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    res = UWORD(0);                                                   \
    for (i = 0; i < (len); i++)                                       \
        res += (expr1) * (expr2);                                     \
    NMOD_RED(res, res, mod);                                          \
} while(0);

// _DOT2_SPLIT   (2 limbs, splitting at 56 bits, 8-unrolling)
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
// 8-unroll: requires  8 * (mod.n - 1)**2 < 2**128
#define _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
do                                                                    \
{                                                                     \
    ulong u0zz = UWORD(0);                                            \
    ulong u1zz = UWORD(0);                                            \
                                                                      \
    for (i = 0; i+7 < (len); )                                        \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
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
    }                                                                 \
    for ( ; i < (len); i++)                                           \
    {                                                                 \
        ulong s0zz, s1zz;                                             \
        umul_ppmm(s1zz, s0zz, (expr1), (expr2));                      \
        add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);               \
    }                                                                 \
                                                                      \
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
    else if (params.method == _DOT2)                                  \
        _NMOD_VEC_DOT2(res, i, len, expr1, expr2, mod)                \
    else if (params.method == _DOT3_ACC)                              \
        _NMOD_VEC_DOT3_ACC(res, i, len, expr1, expr2, mod)            \
    else if (params.method == _DOT3)                                  \
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

/* dot functions: specific algorithms */
ulong _nmod_vec_dot_pow2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot1(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2_half(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3_acc(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);

#if FLINT_BITS == 64
ulong _nmod_vec_dot2_split(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif  // FLINT_BITS == 64

/* general dot functions */

NMOD_VEC_INLINE ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    if (params.method == _DOT1)
        return _nmod_vec_dot1(vec1, vec2, len, mod);

#if FLINT_BITS == 64
    if (params.method == _DOT2_SPLIT)
        return _nmod_vec_dot2_split(vec1, vec2, len, mod, params.pow2_precomp);
#endif // FLINT_BITS == 64

    if (params.method == _DOT2)
        return _nmod_vec_dot2(vec1, vec2, len, mod);

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc(vec1, vec2, len, mod);

    if (params.method == _DOT3)
        return _nmod_vec_dot3(vec1, vec2, len, mod);

    if (params.method == _DOT2_HALF)
        return _nmod_vec_dot2_half(vec1, vec2, len, mod);

    if (params.method == _DOT_POW2)
    {
#if defined(__AVX2__)
        if (mod.n < (UWORD(1) << (FLINT_BITS / 2)))
            return _nmod_vec_dot1(vec1, vec2, len, mod);
        else  // make sure not to use avx 32-bit mul
#endif // defined(__AVX2__)
            return _nmod_vec_dot_pow2(vec1, vec2, len, mod);
    }

    else  // params.method == _DOT0
        return UWORD(0);
}

ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t);

ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, dot_params_t);

/* some IO functions */
#ifdef FLINT_HAVE_FILE
int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod);
#endif

void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod);

#ifdef __cplusplus
}
#endif

#endif

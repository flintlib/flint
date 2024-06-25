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
#include "nmod.h"  // nmod_mul, nmod_fmma

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
void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod);


/* ---- compute dot parameters ---- */

typedef enum
{
    _DOT0 = 0,           /* len == 0 || mod.n == 1 */
    _DOT1 = 1,           /* 1 limb */
#if (FLINT_BITS == 64)
    _DOT2_SPLIT = 2,     /* 2 limbs, modulus < ~2**30.5 (FLINT_BITS == 64 only) */
#endif  // FLINT_BITS == 64
    _DOT2_HALF = 3,      /* 2 limbs, modulus < 2**(FLINT_BITS/2) */
    _DOT2 = 4,           /* 2 limbs */
    _DOT3_ACC = 5,       /* 3 limbs, modulus allowing some accumulation in 2 limbs */
    _DOT3 = 6,           /* 3 limbs */
    _DOT_POW2 = 7,       /* mod.n is a power of 2 */
} dot_method_t;
// if mod.n is a power of 2, we use _DOT_POW2 in all cases
// otherwise, number of limbs of unreduced dot product can be deduced:
// 1 limb  <=>  method <= _DOT1
// 2 limbs <=>  _DOT1 < method <= _DOT2
// 3 limbs <=>  _DOT2 < method

typedef struct
{
    dot_method_t method;
    ulong pow2_precomp;  /* for splitting: (1L << 56) % mod.n */
} dot_params_t;

// for _DOT2_SPLIT
#if (FLINT_BITS == 64)
#   define DOT_SPLIT_BITS 56
#   define DOT_SPLIT_MASK UWORD(72057594037927935) // (1L << DOT_SPLIT_BITS) - 1
#endif // FLINT_BITS == 64

#define _FIXED_LEN_MOD_BOUNDS(fixedlen, onelimb_bnd, twolimb_bnd) \
        if (len == fixedlen)                                      \
        {                                                         \
            if (mod.n <= UWORD(onelimb_bnd))                      \
                return (dot_params_t) {_DOT1, UWORD(0)};          \
            if (mod.n <= UWORD(twolimb_bnd))                      \
                return (dot_params_t) {_DOT2, UWORD(0)};          \
            return (dot_params_t) {_DOT3, UWORD(0)};              \
        }

FLINT_FORCE_INLINE dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod)
{
    if (len == 0 || mod.n == 1)
        return (dot_params_t) {_DOT0, UWORD(0)};
    if ((mod.n & (mod.n - 1)) == 0)
        return (dot_params_t) {_DOT_POW2, UWORD(0)};
    // from here on len >= 1, n > 1 not power of 2

    // short dot products: we use only _DOT1, _DOT2, _DOT3 in that case
    if (len <= 11)
    {
#if FLINT_BITS == 64
        // 64 bits:  k limbs  <=>  n <= ceil(2**(32*k) / sqrt(len))
        _FIXED_LEN_MOD_BOUNDS(11, 1294981365,  5561902608746059656);
        _FIXED_LEN_MOD_BOUNDS(10, 1358187914,  5833372668713515885);
        _FIXED_LEN_MOD_BOUNDS( 9, 1431655766,  6148914691236517206);
        _FIXED_LEN_MOD_BOUNDS( 8, 1518500250,  6521908912666391107);
        _FIXED_LEN_MOD_BOUNDS( 7, 1623345051,  6972213902555716131);
        _FIXED_LEN_MOD_BOUNDS( 6, 1753413057,  7530851732716320753);
        _FIXED_LEN_MOD_BOUNDS( 5, 1920767767,  8249634742471189718);
        _FIXED_LEN_MOD_BOUNDS( 4, 2147483648,  9223372036854775808);
        _FIXED_LEN_MOD_BOUNDS( 3, 2479700525, 10650232656628343402);
        _FIXED_LEN_MOD_BOUNDS( 2, 3037000500, 13043817825332782213);
#else  // FLINT_BITS == 64
        // 32 bits: k limbs  <=>  n <= ceil(2**(16*k) / sqrt(len))
        _FIXED_LEN_MOD_BOUNDS(11, 19760, 1294981365);
        _FIXED_LEN_MOD_BOUNDS(10, 20725, 1358187914);
        _FIXED_LEN_MOD_BOUNDS( 9, 21846, 1431655766);
        _FIXED_LEN_MOD_BOUNDS( 8, 23171, 1518500250);
        _FIXED_LEN_MOD_BOUNDS( 7, 24771, 1623345051);
        _FIXED_LEN_MOD_BOUNDS( 6, 26755, 1753413057);
        _FIXED_LEN_MOD_BOUNDS( 5, 29309, 1920767767);
        _FIXED_LEN_MOD_BOUNDS( 4, 32768, 2147483648);
        _FIXED_LEN_MOD_BOUNDS( 3, 37838, 2479700525);
        _FIXED_LEN_MOD_BOUNDS( 2, 46341, 3037000500);
#endif  // FLINT_BITS == 64
        // remains len == 1
        if (mod.n <= (UWORD(1) << FLINT_BITS / 2))
            return (dot_params_t) {_DOT1, UWORD(0)};
        return (dot_params_t) {_DOT2, UWORD(0)};
    }

    if (mod.n <= UWORD(1) << (FLINT_BITS / 2)) // implies <= 2 limbs
    {
        const ulong t0 = (mod.n - 1) * (mod.n - 1);
        ulong u1, u0;
        umul_ppmm(u1, u0, t0, len);
        if (u1 == 0)  // 1 limb
            return (dot_params_t) {_DOT1, UWORD(0)};

        // u1 != 0 <=> 2 limbs
#if (FLINT_BITS == 64) // _SPLIT: see end of file for these constraints
        if (mod.n <= UWORD(1515531528) && len <= WORD(380368697))
        {
            ulong pow2_precomp;
            NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);
            return (dot_params_t) {_DOT2_SPLIT, pow2_precomp};
        }
#endif
        return (dot_params_t) {_DOT2_HALF, UWORD(0)};
    }
    // from here on, mod.n > 2**(FLINT_BITS / 2)
    // --> unreduced dot cannot fit in 1 limb

    ulong t2, t1, t0, u1, u0;
    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_ssaaaa(t2, t1, t2, t1, UWORD(0), u1);

    if (t2 == 0) // 2 limbs
        return (dot_params_t) {_DOT2, UWORD(0)};

    // 3 limbs:
#if (FLINT_BITS == 64)
    if (mod.n <= UWORD(6521908912666391107))  // room for accumulating 8 terms
#else
    if (mod.n <= UWORD(1518500250))           // room for accumulating 8 terms
#endif
        return (dot_params_t) {_DOT3_ACC, UWORD(0)};

    return (dot_params_t) {_DOT3, UWORD(0)};
}

#undef _FIXED_LEN_MOD_BOUNDS

int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);
int _nmod_vec_dot_bound_limbs_from_params(slong len, nmod_t mod, dot_params_t params);


/* ------ dot product, specific algorithms ------ */

/* vec1[i] * vec2[i] */
ulong _nmod_vec_dot_pow2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot1(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2_half(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3_acc(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
ulong _nmod_vec_dot3(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod);
#if FLINT_BITS == 64
ulong _nmod_vec_dot2_split(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp);
#endif  // FLINT_BITS == 64

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

// * supports 1 <= len <= 11, requires method==DOT1|DOT2|DOT3|DOT_POW2
// * i must be already initialized at the first wanted value
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

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc(vec1, vec2, len, mod);

    if (params.method == _DOT3)
        return _nmod_vec_dot3(vec1, vec2, len, mod);

    if (params.method == _DOT2_HALF)
        return _nmod_vec_dot2_half(vec1, vec2, len, mod);

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

    if (params.method == _DOT2)
        return _nmod_vec_dot2_ptr(vec1, vec2, offset, len, mod);

    if (params.method == _DOT3_ACC)
        return _nmod_vec_dot3_acc_ptr(vec1, vec2, offset, len, mod);

    if (params.method == _DOT3)
        return _nmod_vec_dot3_ptr(vec1, vec2, offset, len, mod);

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


#ifdef __cplusplus
}
#endif

#endif

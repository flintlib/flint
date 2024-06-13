/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

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

int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);


#define DOT_SPLIT_BITS 56
#define DOT_SPLIT_MASK UWORD(72057594037927935) // (1L << DOT_SPLIT_BITS) - 1

// new general dot macro
#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)          \
do                                                                    \
{                                                                     \
    if (nlimbs == 1)                                                  \
    {                                                                 \
        res = UWORD(0);                                               \
        for (i = 0; i < (len); i++)                                   \
            res += (expr1) * (expr2);                                 \
        NMOD_RED(res, res, mod);                                      \
    }                                                                 \
    else if (mod.n <= UWORD(1515531528) && (len) <= WORD(134744072))  \
    {                                                                 \
        ulong dp_lo = 0;                                              \
        unsigned int dp_hi = 0;                                       \
                                                                      \
        for (i = 0; i+7 < (len); )                                    \
        {                                                             \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
            dp_lo += (expr1) * (expr2); i++;                          \
                                                                      \
            dp_hi += dp_lo >> DOT_SPLIT_BITS;                         \
            dp_lo &= DOT_SPLIT_MASK;                                  \
        }                                                             \
                                                                      \
        for ( ; i < (len); i++)                                       \
            dp_lo += (expr1) * (expr2);                               \
                                                                      \
        unsigned int red_pow;                                         \
        NMOD_RED(red_pow, (UWORD(1) << DOT_SPLIT_BITS), mod);         \
        res = (ulong)red_pow * dp_hi + dp_lo;                         \
        NMOD_RED(res, res, mod);                                      \
    }                                                                 \
    else if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))                 \
    {                                                                 \
        ulong s0zz = UWORD(0);                                        \
        ulong s1zz = UWORD(0);                                        \
        for (i = 0; i < (len); i++)                                   \
        {                                                             \
            const ulong prodzz = (expr1) * (expr2);                   \
            add_ssaaaa(s1zz, s0zz, s1zz, s0zz, 0, prodzz);            \
        }                                                             \
        NMOD2_RED2(res, s1zz, s0zz, mod);                             \
    }                                                                 \
    else if (nlimbs == 2)                                             \
    {                                                                 \
        ulong u0zz = UWORD(0);                                        \
        ulong u1zz = UWORD(0);                                        \
                                                                      \
        for (i = 0; i+7 < (len); )                                    \
        {                                                             \
            ulong s0zz, s1zz;                                         \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
            i++;                                                      \
        }                                                             \
        for ( ; i < (len); i++)                                       \
        {                                                             \
            ulong s0zz, s1zz;                                         \
            umul_ppmm(s1zz, s0zz, (expr1), (expr2));                  \
            add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);           \
        }                                                             \
                                                                      \
        NMOD2_RED2(res, u1zz, u0zz, mod);                             \
    }                                                                 \
    else if (nlimbs == 3)                                             \
    {                                                                 \
        ulong t2zz = UWORD(0);                                        \
        ulong t1zz = UWORD(0);                                        \
        ulong t0zz = UWORD(0);                                        \
                                                                      \
        /* we can accumulate 8 terms if n == mod.n is such that */    \
        /*      8 * (n-1)**2 < 2**128, this is equivalent to    */    \
        /*      n <= ceil(sqrt(2**125)) = 6521908912666391107   */    \
        if (mod.n <= 6521908912666391107L)                            \
        {                                                             \
            for (i = 0; i+7 < (len); )                                \
            {                                                         \
                ulong s0zz, s1zz;                                     \
                ulong u0zz = UWORD(0);                                \
                ulong u1zz = UWORD(0);                                \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
                i++;                                                  \
                add_sssaaaaaa(t2zz, t1zz, t0zz,                       \
                              t2zz, t1zz, t0zz,                       \
                              UWORD(0), u1zz, u0zz);                  \
            }                                                         \
                                                                      \
            ulong s0zz, s1zz;                                         \
            ulong u0zz = UWORD(0);                                    \
            ulong u1zz = UWORD(0);                                    \
            for ( ; i < (len); i++)                                   \
            {                                                         \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_ssaaaa(u1zz, u0zz, u1zz, u0zz, s1zz, s0zz);       \
            }                                                         \
                                                                      \
            add_sssaaaaaa(t2zz, t1zz, t0zz,                           \
                          t2zz, t1zz, t0zz,                           \
                          UWORD(0), u1zz, u0zz);                      \
        }                                                             \
        else                                                          \
        {                                                             \
            for (i = 0; i < (len); i++)                               \
            {                                                         \
                ulong s0zz, s1zz;                                     \
                umul_ppmm(s1zz, s0zz, (expr1), (expr2));              \
                add_sssaaaaaa(t2zz, t1zz, t0zz,                       \
                              t2zz, t1zz, t0zz,                       \
                              UWORD(0), s1zz, s0zz);                  \
            }                                                         \
        }                                                             \
                                                                      \
        NMOD_RED(t2zz, t2zz, mod);                                    \
        NMOD_RED3(res, t2zz, t1zz, t0zz, mod);                        \
    }                                                                 \
    else   /* nlimbs == 0 */                                          \
    {                                                                 \
        res = UWORD(0);                                               \
    }                                                                 \
} while(0);

ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs);
ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs);

ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, int nlimbs);

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

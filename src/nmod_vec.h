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

void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec,
                                        slong len, nmod_t mod);

void _nmod_vec_add(nn_ptr res, nn_srcptr vec1,
                        nn_srcptr vec2, slong len, nmod_t mod);

void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1,
                        nn_srcptr vec2, slong len, nmod_t mod);

void _nmod_vec_neg(nn_ptr res, nn_srcptr vec,
                                            slong len, nmod_t mod);

void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec,
                            slong len, ulong c, nmod_t mod);

void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec,
                            slong len, ulong c, nmod_t mod);

void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec,
                            slong len, ulong c, nmod_t mod);

int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);


#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)                \
    do                                                                      \
    {                                                                       \
        ulong s0, s1, s2, t0, t1;                                       \
        s0 = s1 = s2 = UWORD(0);                                            \
        switch (nlimbs)                                                     \
        {                                                                   \
            case 1:                                                         \
                for (i = 0; i < (len); i++)                                 \
                {                                                           \
                    s0 += (expr1) * (expr2);                                \
                }                                                           \
                NMOD_RED(s0, s0, mod);                                      \
                break;                                                      \
            case 2:                                                         \
                if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))                \
                {                                                           \
                    for (i = 0; i < (len); i++)                             \
                    {                                                       \
                        t0 = (expr1) * (expr2);                             \
                        add_ssaaaa(s1, s0, s1, s0, 0, t0);                  \
                    }                                                       \
                }                                                           \
                else if ((len) < 8)                                         \
                {                                                           \
                    for (i = 0; i < len; i++)                               \
                    {                                                       \
                        umul_ppmm(t1, t0, (expr1), (expr2));                \
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);                 \
                    }                                                       \
                }                                                           \
                else                                                        \
                {                                                           \
                    ulong v0, v1, u0, u1;                               \
                    i = 0;                                                  \
                    if ((len) & 1)                                          \
                        umul_ppmm(v1, v0, (expr1), (expr2));                \
                    else                                                    \
                        v0 = v1 = 0;                                        \
                    for (i = (len) & 1; i < (len); i++)                     \
                    {                                                       \
                        umul_ppmm(t1, t0, (expr1), (expr2));                \
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);                 \
                        i++;                                                \
                        umul_ppmm(u1, u0, (expr1), (expr2));                \
                        add_ssaaaa(v1, v0, v1, v0, u1, u0);                 \
                    }                                                       \
                    add_ssaaaa(s1, s0, s1, s0, v1, v0);                     \
                }                                                           \
                NMOD2_RED2(s0, s1, s0, mod);                                \
                break;                                                      \
            default:                                                        \
                for (i = 0; i < (len); i++)                                 \
                {                                                           \
                    umul_ppmm(t1, t0, (expr1), (expr2));                    \
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);       \
                }                                                           \
                NMOD_RED(s2, s2, mod);                                      \
                NMOD_RED3(s0, s2, s1, s0, mod);                             \
                break;                                                      \
        }                                                                   \
        res = s0;                                                           \
    } while (0);

ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset,
    slong len, nmod_t mod, int nlimbs);

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

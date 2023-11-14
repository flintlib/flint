/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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
mp_ptr _nmod_vec_init(slong len)
{
   return (mp_ptr) flint_malloc(len * sizeof(mp_limb_t));
}

NMOD_VEC_INLINE
void _nmod_vec_clear(mp_ptr vec)
{
   flint_free(vec);
}

void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, slong len, nmod_t mod);

NMOD_VEC_INLINE
void _nmod_vec_zero(mp_ptr vec, slong len)
{
   flint_mpn_zero(vec, len);
}

flint_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len);

NMOD_VEC_INLINE
void _nmod_vec_set(mp_ptr res, mp_srcptr vec, slong len)
{
   flint_mpn_copyi(res, vec, len);
}

NMOD_VEC_INLINE
void _nmod_vec_swap(mp_ptr a, mp_ptr b, slong length)
{
    slong i;
    for (i = 0; i < length; i++)
    {
        mp_limb_t t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}

NMOD_VEC_INLINE
int _nmod_vec_equal(mp_srcptr vec, mp_srcptr vec2, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != vec2[i]) return 0;

   return 1;
}

NMOD_VEC_INLINE
int _nmod_vec_is_zero(mp_srcptr vec, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != 0) return 0;

   return 1;
}

void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec,
                                        slong len, nmod_t mod);

void _nmod_vec_add(mp_ptr res, mp_srcptr vec1,
                        mp_srcptr vec2, slong len, nmod_t mod);

void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1,
                        mp_srcptr vec2, slong len, nmod_t mod);

void _nmod_vec_neg(mp_ptr res, mp_srcptr vec,
                                            slong len, nmod_t mod);

void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec,
                            slong len, mp_limb_t c, nmod_t mod);

void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec,
                            slong len, mp_limb_t c, nmod_t mod);

void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec,
                            slong len, mp_limb_t c, nmod_t mod);

int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);


#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)                \
    do                                                                      \
    {                                                                       \
        mp_limb_t s0, s1, s2, t0, t1;                                       \
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
                    mp_limb_t v0, v1, u0, u1;                               \
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

mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

mp_limb_t _nmod_vec_dot_rev(mp_srcptr vec1, mp_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset,
    slong len, nmod_t mod, int nlimbs);

/* some IO functions */
#ifdef FLINT_HAVE_FILE
int _nmod_vec_fprint_pretty(FILE * file, mp_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_fprint(FILE * f, mp_srcptr vec, slong len, nmod_t mod);
#endif

void _nmod_vec_print_pretty(mp_srcptr vec, slong len, nmod_t mod);
int _nmod_vec_print(mp_srcptr vec, slong len, nmod_t mod);



#ifdef __cplusplus
}
#endif

#endif


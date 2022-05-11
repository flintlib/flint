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

#include "nmod_mini.h"
#include "nmod_vec_mini.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NMOD_VEC_NORM(vec, i)   \
    do                          \
    {                           \
        while ((i) && vec[(i) - 1] == UWORD(0)) \
            (i)--;              \
    } while (0)

NMOD_VEC_INLINE
ulong_ptr _nmod_vec_init(slong len)
{
   return flint_malloc(len * sizeof(ulong));
}

NMOD_VEC_INLINE
void _nmod_vec_clear(ulong_ptr vec)
{
   flint_free(vec);
}

FLINT_DLL void _nmod_vec_randtest(ulong_ptr vec, flint_rand_t state, slong len, nmod_t mod);

#define _NMOD_VEC_SET   FLINT_MPN_COPYI
#define _NMOD_VEC_ZERO  FLINT_MPN_ZERO

/* FIXME: Write these two functions. These are currently commented as the user
 * should not be able to access flint-impl.h. */
/* NMOD_VEC_INLINE */
/* void _nmod_vec_set(ulong_ptr res, ulong_srcptr vec, slong len) */
/* { */
/*    FLINT_MPN_COPYI(res, vec, len); */
/* } */

/* NMOD_VEC_INLINE */
/* void _nmod_vec_zero(ulong_ptr vec, slong len) */
/* { */
/*    FLINT_MPN_ZERO(vec, len); */
/* } */

NMOD_VEC_INLINE
void _nmod_vec_swap(ulong_ptr a, ulong_ptr b, slong length)
{
    slong i;
    for (i = 0; i < length; i++)
    {
        ulong t = a[i];
        a[i] = b[i];
        b[i] = t;
    }
}

FLINT_DLL void _nmod_vec_reduce(ulong_ptr res, ulong_srcptr vec, 
                                        slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_add(ulong_ptr res, ulong_srcptr vec1, 
                        ulong_srcptr vec2, slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_sub(ulong_ptr res, ulong_srcptr vec1, 
                        ulong_srcptr vec2, slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_neg(ulong_ptr res, ulong_srcptr vec, 
                                            slong len, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_mul_nmod(ulong_ptr res, ulong_srcptr vec, 
                            slong len, ulong c, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_mul_nmod_shoup(ulong_ptr res, ulong_srcptr vec, 
                            slong len, ulong c, nmod_t mod);

FLINT_DLL void _nmod_vec_scalar_addmul_nmod(ulong_ptr res, ulong_srcptr vec, 
                            slong len, ulong c, nmod_t mod);

FLINT_DLL int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod);


#define NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)    \
    do                                                          \
    {                                                           \
        ulong s0, s1, s2, t0, t1;                               \
        s0 = s1 = s2 = UWORD(0);                                \
        switch (nlimbs)                                         \
        {                                                       \
            case 1:                                             \
                for (i = 0; i < (len); i++)                     \
                {                                               \
                    s0 += (expr1) * (expr2);                    \
                }                                               \
                NMOD_RED(s0, s0, mod);                          \
                break;                                          \
            case 2:                                             \
                if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))    \
                {                                               \
                    for (i = 0; i < (len); i++)                 \
                    {                                           \
                        t0 = (expr1) * (expr2);                 \
                        add_ssaaaa(s1, s0, s1, s0, 0, t0);      \
                    }                                           \
                }                                               \
                else if ((len) < 8)                             \
                {                                               \
                    for (i = 0; i < len; i++)                   \
                    {                                           \
                        umul_ppmm(t1, t0, (expr1), (expr2));    \
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);     \
                    }                                           \
                }                                               \
                else                                            \
                {                                               \
                    ulong v0, v1, u0, u1;                       \
                    i = 0;                                      \
                    if ((len) & 1)                              \
                        umul_ppmm(v1, v0, (expr1), (expr2));    \
                    else                                        \
                        v0 = v1 = 0;                            \
                    for (i = (len) & 1; i < (len); i++)         \
                    {                                           \
                        umul_ppmm(t1, t0, (expr1), (expr2));    \
                        add_ssaaaa(s1, s0, s1, s0, t1, t0);     \
                        i++;                                    \
                        umul_ppmm(u1, u0, (expr1), (expr2));    \
                        add_ssaaaa(v1, v0, v1, v0, u1, u0);     \
                    }                                           \
                    add_ssaaaa(s1, s0, s1, s0, v1, v0);         \
                }                                               \
                NMOD2_RED2(s0, s1, s0, mod);                    \
                break;                                          \
            default:                                            \
                for (i = 0; i < (len); i++)                     \
                {                                               \
                    umul_ppmm(t1, t0, (expr1), (expr2));        \
                    add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0); \
                }                                               \
                NMOD_RED(s2, s2, mod);                          \
                NMOD_RED3(s0, s2, s1, s0, mod);                 \
                break;                                          \
        }                                                       \
        res = s0;                                               \
    } while (0);

FLINT_DLL ulong _nmod_vec_dot(ulong_srcptr vec1, ulong_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

FLINT_DLL ulong _nmod_vec_dot_rev(ulong_srcptr vec1, ulong_srcptr vec2,
    slong len, nmod_t mod, int nlimbs);

FLINT_DLL ulong _nmod_vec_dot_ptr(ulong_srcptr vec1, const ulong_ptr * vec2, slong offset,
    slong len, nmod_t mod, int nlimbs);

#ifdef __cplusplus
}
#endif

#endif

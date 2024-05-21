/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_mod_vec.h"

#if defined(__GNUC__)
# pragma GCC push_options
# pragma GCC optimize ("unroll-loops")
#endif

#if TUNE_PROGRAM && N_MOD_VEC_ADD_METHOD == 0
/* Corresponds to _nmod_add */
# define N_MOD_VEC_ADD(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    for (ix = 0; ix < len; ix++)                    \
    {                                               \
        ulong sum = up[ix] + vp[ix];                \
        rp[ix] = sum - mod + ((((slong) (sum - mod)) >> (FLINT_BITS - 1)) & mod); \
    }                                               \
} while (0)
#elif N_MOD_VEC_ADD_METHOD == 0
/* Corresponds to _nmod_add and nmod_add */
# define N_MOD_VEC_ADD(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    if ((mod >> (FLINT_BITS - 1)) == 0)             \
        for (ix = 0; ix < len; ix++)                \
        {                                           \
            ulong sum = up[ix] + vp[ix];            \
            rp[ix] = sum - mod + ((((slong) (sum - mod)) >> (FLINT_BITS - 1)) & mod); \
        }                                           \
    else                                            \
        for (ix = 0; ix < len; ix++)                \
        {                                           \
            ulong neg = mod - up[ix];               \
            rp[ix] = (neg > vp[ix]) ? up[ix] + vp[ix] : vp[ix] - neg; \
        }                                           \
} while (0)
#elif N_MOD_VEC_ADD_METHOD == 1
/* Corresponds to nmod_add */
# define N_MOD_VEC_ADD(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    for (ix = 0; ix < len; ix++)                    \
    {                                               \
        ulong neg = mod - up[ix];                   \
        rp[ix] = (neg > vp[ix]) ? up[ix] + vp[ix] : vp[ix] - neg; \
    }                                               \
} while (0)
#else
# error
#endif

#if TUNE_PROGRAM && N_MOD_VEC_SUB_METHOD == 0
/* Corresponds to _nmod_sub */
# define N_MOD_VEC_SUB(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    for (ix = 0; ix < len; ix++)                    \
    {                                               \
        ulong diff = up[ix] - vp[ix];               \
        rp[ix] = diff + ((((slong) diff) >> (FLINT_BITS - 1)) & mod); \
    }                                               \
} while (0)
#elif N_MOD_VEC_SUB_METHOD == 0
/* Corresponds to _nmod_sub and nmod_sub */
# define N_MOD_VEC_SUB(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    if ((mod >> (FLINT_BITS - 1)) == 0)             \
        for (ix = 0; ix < len; ix++)                \
        {                                           \
            ulong diff = up[ix] - vp[ix];           \
            rp[ix] = diff + ((((slong) diff) >> (FLINT_BITS - 1)) & mod); \
        }                                           \
    else                                            \
        for (ix = 0; ix < len; ix++)                \
        {                                           \
            ulong diff = up[ix] - vp[ix];           \
            rp[ix] = (up[ix] < vp[ix]) ? mod + diff : diff; \
        }                                           \
} while (0)
#elif N_MOD_VEC_SUB_METHOD == 1
/* Corresponds to nmod_sub */
# define N_MOD_VEC_SUB(rp, up, vp, len, mod)        \
do {                                                \
    slong ix;                                       \
    for (ix = 0; ix < len; ix++)                    \
    {                                               \
        ulong diff = up[ix] - vp[ix];               \
        rp[ix] = (up[ix] < vp[ix]) ? mod + diff : diff; \
    }                                               \
} while (0)
#else
# error
#endif

void
_n_mod_vec_add(nn_ptr restrict rp, nn_srcptr restrict up, nn_srcptr restrict vp, slong len, ulong mod)
{
    if (len <= 0)
        FLINT_UNREACHABLE;

    N_MOD_VEC_ADD(rp, up, vp, len, mod);
}

void
_n_mod_vec_sub(nn_ptr restrict rp, nn_srcptr restrict up, nn_srcptr restrict vp, slong len, ulong mod)
{
    if (len <= 0)
        FLINT_UNREACHABLE;

    N_MOD_VEC_SUB(rp, up, vp, len, mod);
}

#if defined(__GNUC__)
# pragma GCC pop_options
#endif

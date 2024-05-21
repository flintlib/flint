/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef N_MOD_VEC_H
#define N_MOD_VEC_H

#ifdef N_MOD_VEC_INLINES_C
# define N_MOD_VEC_INLINE
#else
# define N_MOD_VEC_INLINE static inline
#endif

#include "n_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* memory management *********************************************************/

N_MOD_VEC_INLINE nn_ptr _n_mod_vec_init(slong len)
{
    return (nn_ptr) flint_malloc(len * sizeof(ulong));
}

N_MOD_VEC_INLINE void _n_mod_vec_clear(nn_ptr up)
{
    flint_free(up);
}

void _n_mod_vec_is_canonical(nn_srcptr, slong, n_mod_ctx_srcptr);

/* randomisation *************************************************************/

void _n_mod_vec_rand(nn_ptr, flint_rand_t, slong, n_mod_ctx_srcptr);

/* assignments ***************************************************************/

N_MOD_VEC_INLINE void _n_mod_vec_zero(nn_ptr up, slong len)
{
#if defined(__GNUC__)
    __builtin_memset(up, 0, sizeof(ulong) * len);
#else
    slong ix;
    for (ix = 0; ix < len; ix++)
        up[ix] = 0;
#endif
}

N_MOD_VEC_INLINE
void _n_mod_vec_set(nn_ptr rp, nn_srcptr ip, slong len)
{
#if defined(__GNUC__)
    __builtin_memcpy(rp, ip, sizeof(ulong) * len);
#else
    slong ix;
    for (ix = 0; ix < len; ix++)
        rp[ix] = ip[ix];
#endif
}

/* comparisons ***************************************************************/

N_MOD_VEC_INLINE
int _n_mod_vec_equal(nn_srcptr up, nn_srcptr vp, slong len)
{
#if defined(__GNUC__)
    return __builtin_memcmp(up, vp, sizeof(ulong) * len) != 0;
#else
    slong ix;

    for (ix = 0; ix < len; ix++)
        if (up[ix] != vp[ix])
            return 0;

    return 1;
#endif
}

N_MOD_VEC_INLINE
int _n_mod_vec_is_zero(nn_srcptr up, slong len)
{
    slong ix;

    for (ix = 0; ix < len; ix++)
        if (up[ix] != UWORD(0))
            return 0;

    return 1;
}

/* addition ******************************************************************/

void _n_mod_vec_neg(nn_ptr, nn_srcptr, slong, ulong);

void _n_mod_vec_add(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);
void _n_mod_vec_sub(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);

/* dot products **************************************************************/

/*
   0: unreduced dot product fits inside a word
   1: mod < 2^{FLINT_BITS / 2} / 16
   2: mod < 2^{FLINT_BITS / 2}
   3: mod < 2^{FLINT_BITS} / 16
   4: mod < 2^{FLINT_BITS}
*/
ulong _n_mod_vec_dot_0(nn_srcptr, nn_srcptr, slong, n_mod_ctx_srcptr);
ulong _n_mod_vec_dot_1(nn_srcptr, nn_srcptr, slong, n_mod_ctx_srcptr);
ulong _n_mod_vec_dot_2(nn_srcptr, nn_srcptr, slong, n_mod_ctx_srcptr);
ulong _n_mod_vec_dot_3(nn_srcptr, nn_srcptr, slong, n_mod_ctx_srcptr);
ulong _n_mod_vec_dot_4(nn_srcptr, nn_srcptr, slong, n_mod_ctx_srcptr);

#ifdef __cplusplus
}
#endif

#endif

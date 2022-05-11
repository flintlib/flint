/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_VEC_MINI_H
#define NMOD_VEC_MINI_H

#ifdef NMOD_VEC_INLINES_C
#define NMOD_VEC_INLINE FLINT_DLL
#else
#define NMOD_VEC_INLINE static __inline__
#endif

#include "nmod_mini.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL flint_bitcnt_t _nmod_vec_max_bits(ulong_srcptr vec, slong len);

NMOD_VEC_INLINE
int _nmod_vec_equal(ulong_srcptr vec, ulong_srcptr vec2, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != vec2[i]) return 0;

   return 1;
}

NMOD_VEC_INLINE
int _nmod_vec_is_zero(ulong_srcptr vec, slong len)
{
   slong i;

   for (i = 0; i < len; i++)
      if (vec[i] != 0) return 0;

   return 1;
}

#ifdef __cplusplus
}
#endif

#endif

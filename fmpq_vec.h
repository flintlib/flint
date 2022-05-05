/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPQ_VEC_H
#define FMPQ_VEC_H

#ifdef FMPQ_VEC_INLINES_C
#define FMPQ_VEC_INLINE FLINT_DLL
#else
#define FMPQ_VEC_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

FLINT_DLL fmpq * _fmpq_vec_init(slong len);

FLINT_DLL void _fmpq_vec_clear(fmpq * vec, slong len);

/*  Randomisation  ***********************************************************/

FLINT_DLL void _fmpq_vec_randtest(fmpq * f, flint_rand_t state, 
                        slong len, flint_bitcnt_t bits);

FLINT_DLL void _fmpq_vec_randtest_uniq_sorted(fmpq * vec,
                        flint_rand_t state, slong len, flint_bitcnt_t bits);

/* Sorting  ******************************************************************/

FLINT_DLL void _fmpq_vec_sort(fmpq * vec, slong len);

/*  Conversions  *************************************************************/

FLINT_DLL void _fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len);

/*  Dot product  **************************************************/

FLINT_DLL void _fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len);

/*  Input and output  ********************************************************/

#if defined (FILE)                  \
  || defined (H_STDIO)              \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)               \
  || defined (__DEFINED_FILE)
FLINT_DLL int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len);

FMPQ_VEC_INLINE
int _fmpq_vec_print(const fmpq * vec, slong len)
{
    return _fmpq_vec_fprint(stdout, vec, len);
}
#endif

#ifdef __cplusplus
}
#endif

#endif

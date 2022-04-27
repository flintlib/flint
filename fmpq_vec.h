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

#include <gmp.h>
#include "fmpq.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

/* _fmpq_vec_init and _fmpq_vec_clear are declared in fmpq.h for backward
   compatibility */

/*  Randomisation  ***********************************************************/

FLINT_DLL void _fmpq_vec_randtest(fmpq * f, flint_rand_t state, 
                        slong len, flint_bitcnt_t bits);

FLINT_DLL void _fmpq_vec_randtest_uniq_sorted(fmpq * vec,
                        flint_rand_t state, slong len, flint_bitcnt_t bits);

/* Sorting  ******************************************************************/

FLINT_DLL void _fmpq_vec_sort(fmpq * vec, slong len);

/*  Conversions  *************************************************************/

FLINT_DLL void _fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len);

FLINT_DLL void _fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den,
                                                    const fmpq * a, slong len);

FLINT_DLL void _fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den,
                                                    const fmpq * a, slong len);

/*  Dot product  **************************************************/

FLINT_DLL void _fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len);

/*  Input and output  ********************************************************/

FLINT_DLL int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len);

FMPQ_VEC_INLINE
int _fmpq_vec_print(const fmpq * vec, slong len)
{
    return _fmpq_vec_fprint(stdout, vec, len);
}


#ifdef __cplusplus
}
#endif

#endif

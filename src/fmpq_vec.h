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
#define FMPQ_VEC_INLINE
#else
#define FMPQ_VEC_INLINE static inline
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Memory management  *******************************************************/

fmpq * _fmpq_vec_init(slong len);
void _fmpq_vec_clear(fmpq * vec, slong len);

/*  Randomisation  ***********************************************************/

void _fmpq_vec_randtest(fmpq * f, flint_rand_t state,
                        slong len, flint_bitcnt_t bits);

void _fmpq_vec_randtest_uniq_sorted(fmpq * vec,
                        flint_rand_t state, slong len, flint_bitcnt_t bits);

/* Sorting  ******************************************************************/

void _fmpq_vec_sort(fmpq * vec, slong len);

/*  Conversions  *************************************************************/

void _fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len);

void _fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den,
                                                    const fmpq * a, slong len);

void _fmpq_vec_get_fmpz_vec_fmpz(fmpz* num, fmpz_t den,
                                                    const fmpq * a, slong len);

/*  Dot product  **************************************************/

void _fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len);

/*  Input and output  ********************************************************/

#ifdef FLINT_HAVE_FILE
int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len);
#endif

int _fmpq_vec_print(const fmpq * vec, slong len);

#ifdef __cplusplus
}
#endif

#endif

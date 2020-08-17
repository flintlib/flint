/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPFR_MAT_H
#define MPFR_MAT_H

#ifdef MPFR_MAT_INLINES_C
#define MPFR_MAT_INLINE FLINT_DLL
#else
#define MPFR_MAT_INLINE static __inline__
#endif

#include <gmp.h>
#include <mpfr.h> 

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    __mpfr_struct * entries;
    slong r;
    slong c;
    flint_bitcnt_t prec;
    __mpfr_struct ** rows;
} mpfr_mat_struct;

/* fmpz_mat_t allows reference-like semantics for fmpz_mat_struct */
typedef mpfr_mat_struct mpfr_mat_t[1];

MPFR_MAT_INLINE
__mpfr_struct * mpfr_mat_entry(const mpfr_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FLINT_DLL void mpfr_mat_init(mpfr_mat_t mat, slong rows, slong cols, mpfr_prec_t prec);

FLINT_DLL void mpfr_mat_swap(mpfr_mat_t mat1, mpfr_mat_t mat2);

FLINT_DLL void mpfr_mat_set(mpfr_mat_t mat1, const mpfr_mat_t mat2);

FLINT_DLL void mpfr_mat_clear(mpfr_mat_t mat);

FLINT_DLL int mpfr_mat_equal(const mpfr_mat_t mat1, const mpfr_mat_t mat2);

FLINT_DLL void mpfr_mat_zero(mpfr_mat_t mat);

/* Random matrix generation  *************************************************/

FLINT_DLL void mpfr_mat_randtest(mpfr_mat_t mat, flint_rand_t state);

/* Multiplication */

FLINT_DLL void mpfr_mat_mul_classical(mpfr_mat_t C, const mpfr_mat_t A, const mpfr_mat_t B, mpfr_rnd_t rnd);

#ifdef __cplusplus
}
#endif

#endif


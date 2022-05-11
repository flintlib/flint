/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MAT_MINI_H
#define FMPZ_MAT_MINI_H

#ifdef FMPZ_MAT_INLINES_C
#define FMPZ_MAT_INLINE FLINT_DLL
#else
#define FMPZ_MAT_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* element access *************************************************************/

FMPZ_MAT_INLINE
fmpz * fmpz_mat_entry(const fmpz_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FMPZ_MAT_INLINE
slong fmpz_mat_nrows(const fmpz_mat_t mat)
{
   return mat->r;
}

FMPZ_MAT_INLINE
slong fmpz_mat_ncols(const fmpz_mat_t mat)
{
   return mat->c;
}

/* assignment *****************************************************************/

FLINT_DLL void fmpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2);

FLINT_DLL void fmpz_mat_zero(fmpz_mat_t mat);
FLINT_DLL void fmpz_mat_one(fmpz_mat_t mat);

/* comparison *****************************************************************/

FLINT_DLL int fmpz_mat_equal(const fmpz_mat_t mat1, const fmpz_mat_t mat2);

FLINT_DLL int fmpz_mat_is_zero(const fmpz_mat_t mat);
FLINT_DLL int fmpz_mat_is_one(const fmpz_mat_t mat);

FMPZ_MAT_INLINE
int fmpz_mat_is_empty(const fmpz_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

FMPZ_MAT_INLINE
int fmpz_mat_is_square(const fmpz_mat_t mat)
{
    return (mat->r == mat->c);
}

/* random generation **********************************************************/

FLINT_DLL void fmpz_mat_randops(fmpz_mat_t mat, flint_rand_t state, slong count);

/* concatenation **************************************************************/

FLINT_DLL void fmpz_mat_concat_horizontal(fmpz_mat_t res, const fmpz_mat_t mat1,  const fmpz_mat_t mat2);
FLINT_DLL void fmpz_mat_concat_vertical(fmpz_mat_t res, const fmpz_mat_t mat1,  const fmpz_mat_t mat2);

#ifdef __cplusplus
}
#endif

#endif

/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MPF_MAT_H
#define MPF_MAT_H

#ifdef MPF_MAT_INLINES_C
#define MPF_MAT_INLINE
#else
#define MPF_MAT_INLINE static inline
#endif

#include "fmpz_types.h"
#include "mpf_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    mpf * entries;
    slong r;
    slong c;
    flint_bitcnt_t prec;
    mpf ** rows;
} mpf_mat_struct;

typedef mpf_mat_struct mpf_mat_t[1];

MPF_MAT_INLINE
mpf * mpf_mat_entry(const mpf_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

MPF_MAT_INLINE
slong mpf_mat_nrows(const mpf_mat_t mat)
{
    return mat->r;
}

MPF_MAT_INLINE
slong mpf_mat_ncols(const mpf_mat_t mat)
{
    return mat->c;
}

/* Memory management  ********************************************************/

void mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, flint_bitcnt_t prec);

void mpf_mat_swap(mpf_mat_t mat1, mpf_mat_t mat2);

MPF_MAT_INLINE void
mpf_mat_swap_entrywise(mpf_mat_t mat1, mpf_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < mpf_mat_nrows(mat1); i++)
        for (j = 0; j < mpf_mat_ncols(mat1); j++)
            mpf_swap(mpf_mat_entry(mat2, i, j), mpf_mat_entry(mat1, i, j));
}

void mpf_mat_set(mpf_mat_t mat1, const mpf_mat_t mat2);

void mpf_mat_clear(mpf_mat_t mat);

int mpf_mat_equal(const mpf_mat_t mat1, const mpf_mat_t mat2);

int mpf_mat_approx_equal(const mpf_mat_t mat1, const mpf_mat_t mat2, flint_bitcnt_t bits);

int mpf_mat_is_zero(const mpf_mat_t mat);

MPF_MAT_INLINE int
mpf_mat_is_empty(const mpf_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

MPF_MAT_INLINE int
mpf_mat_is_square(const mpf_mat_t mat)
{
    return (mat->r == mat->c);
}

void mpf_mat_zero(mpf_mat_t mat);

void mpf_mat_one(mpf_mat_t mat);


/* Conversions ***************************************************************/

void mpf_mat_set_fmpz_mat(mpf_mat_t B, const fmpz_mat_t A);

/* Backwards compatibility (this will generate errors for non-GCC compatible
 * compilers) */
#define fmpz_mat_get_mpf_mat(B, A) _Pragma("GCC warning \"'fmpz_mat_get_mpf_mat' is deprecated in favor for 'mpf_mat_set_fmpz_mat'\"") mpf_mat_set_fmpz_mat(B, A)

/* Input and output  *********************************************************/

void mpf_mat_print(const mpf_mat_t mat);

/* Random matrix generation  *************************************************/

void mpf_mat_randtest(mpf_mat_t mat, flint_rand_t state, flint_bitcnt_t bits);

/* Multiplication */

void mpf_mat_mul(mpf_mat_t C, const mpf_mat_t A, const mpf_mat_t B);

/* Permutations */

MPF_MAT_INLINE void
mpf_mat_swap_rows(mpf_mat_t mat, slong r, slong s)
{
    if (r != s)
    {
        mpf * u;

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

/* Gram-Schmidt Orthogonalisation and QR Decomposition  ********************************************************/

void mpf_mat_gso(mpf_mat_t B, const mpf_mat_t A);

void mpf_mat_qr(mpf_mat_t Q, mpf_mat_t R, const mpf_mat_t A);

#ifdef __cplusplus
}
#endif

#endif


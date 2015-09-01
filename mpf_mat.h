/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#ifndef MPF_MAT_H
#define MPF_MAT_H

#ifdef MPF_MAT_INLINES_C
#define MPF_MAT_INLINE FLINT_DLL
#else
#define MPF_MAT_INLINE static __inline__
#endif

#include <math.h>
#include <stdio.h>
#include "mpf_vec.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    mpf * entries;
    slong r;
    slong c;
    mp_bitcnt_t prec;
    mpf ** rows;
} mpf_mat_struct;

typedef mpf_mat_struct mpf_mat_t[1];

MPF_MAT_INLINE
mpf * mpf_mat_entry(const mpf_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

/* Memory management  ********************************************************/

FLINT_DLL void mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, mp_bitcnt_t prec);

FLINT_DLL void mpf_mat_swap(mpf_mat_t mat1, mpf_mat_t mat2);

FLINT_DLL void mpf_mat_set(mpf_mat_t mat1, const mpf_mat_t mat2);

FLINT_DLL void mpf_mat_clear(mpf_mat_t mat);

FLINT_DLL int mpf_mat_equal(const mpf_mat_t mat1, const mpf_mat_t mat2);

FLINT_DLL int mpf_mat_approx_equal(const mpf_mat_t mat1, const mpf_mat_t mat2, mp_bitcnt_t bits);

FLINT_DLL int mpf_mat_is_zero(const mpf_mat_t mat);

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

FLINT_DLL void mpf_mat_zero(mpf_mat_t mat);

FLINT_DLL void mpf_mat_one(mpf_mat_t mat);

/* Input and output  *********************************************************/

FLINT_DLL void mpf_mat_print(const mpf_mat_t mat);

/* Random matrix generation  *************************************************/

FLINT_DLL void mpf_mat_randtest(mpf_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

/* Multiplication */

FLINT_DLL void mpf_mat_mul(mpf_mat_t C, const mpf_mat_t A, const mpf_mat_t B);

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

FLINT_DLL void mpf_mat_gso(mpf_mat_t B, const mpf_mat_t A);

FLINT_DLL void mpf_mat_qr(mpf_mat_t Q, mpf_mat_t R, const mpf_mat_t A);

#ifdef __cplusplus
}
#endif

#endif


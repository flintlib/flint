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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPQ_MAT_H
#define FMPQ_MAT_H

#ifdef FMPQ_MAT_INLINES_C
#define FMPQ_MAT_INLINE FLINT_DLL
#else
#define FMPQ_MAT_INLINE static __inline__
#endif

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    fmpq * entries;
    slong r;
    slong c;
    fmpq ** rows;
} fmpq_mat_struct;

typedef fmpq_mat_struct fmpq_mat_t[1];

FMPQ_MAT_INLINE
fmpq * fmpq_mat_entry(const fmpq_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FMPQ_MAT_INLINE
fmpz * fmpq_mat_entry_num(const fmpq_mat_t mat, slong i, slong j)
{
   return (fmpz *)(&((*fmpq_mat_entry(mat, i, j)).num));
}

FMPQ_MAT_INLINE
fmpz * fmpq_mat_entry_den(const fmpq_mat_t mat, slong i, slong j)
{
   return (fmpz *)(&((*fmpq_mat_entry(mat, i, j)).den));
}

FMPQ_MAT_INLINE
slong fmpq_mat_nrows(const fmpq_mat_t mat)
{
   return mat->r;
}

FMPQ_MAT_INLINE
slong fmpq_mat_ncols(const fmpq_mat_t mat)
{
   return mat->c;
}

FLINT_DLL void fmpq_mat_init(fmpq_mat_t mat, slong rows, slong cols);

FLINT_DLL void fmpq_mat_clear(fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_swap(fmpq_mat_t mat1, fmpq_mat_t mat2);

/* Windows and concatenation */

FLINT_DLL void fmpq_mat_window_init(fmpq_mat_t window, const fmpq_mat_t mat, slong r1,
    slong c1, slong r2, slong c2);

FLINT_DLL void fmpq_mat_window_clear(fmpq_mat_t window);

FLINT_DLL void fmpq_mat_concat_horizontal(fmpq_mat_t res,
                           const fmpq_mat_t mat1,  const fmpq_mat_t mat2);

FLINT_DLL void fmpq_mat_concat_vertical(fmpq_mat_t res,
                           const fmpq_mat_t mat1,  const fmpq_mat_t mat2);

/* Input and output  */

FLINT_DLL void fmpq_mat_print(const fmpq_mat_t mat);

/* Random matrix generation **************************************************/

FLINT_DLL void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

FLINT_DLL void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

/* Special matrices **********************************************************/

FLINT_DLL void fmpq_mat_hilbert_matrix(fmpq_mat_t mat);

/* Basic assignment **********************************************************/

FLINT_DLL void fmpq_mat_set(fmpq_mat_t dest, const fmpq_mat_t src);

FLINT_DLL void fmpq_mat_zero(fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_one(fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_transpose(fmpq_mat_t rop, const fmpq_mat_t op);

/* Addition, scalar multiplication  ******************************************/

FLINT_DLL void fmpq_mat_add(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2);

FLINT_DLL void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2);

FLINT_DLL void fmpq_mat_neg(fmpq_mat_t rop, const fmpq_mat_t op);

FLINT_DLL void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x);

FLINT_DLL void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x);

/* Basic comparison and properties *******************************************/

FLINT_DLL int fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2);

FLINT_DLL int fmpq_mat_is_integral(const fmpq_mat_t mat);

FLINT_DLL int fmpq_mat_is_zero(const fmpq_mat_t mat);

FMPQ_MAT_INLINE
int fmpq_mat_is_empty(const fmpq_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

FMPQ_MAT_INLINE
int fmpq_mat_is_square(const fmpq_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Integer matrix conversion *************************************************/

FLINT_DLL int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, const fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den,
    const fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den,
    const fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

FLINT_DLL void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2,
        fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2);

FLINT_DLL void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, const fmpq_mat_t mat,
    const fmpz_t mod);

FLINT_DLL void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src);

FLINT_DLL void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod,
    const fmpz_t div);

FLINT_DLL int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod,
    const fmpz_t mod);

/* Matrix multiplication *****************************************************/

FLINT_DLL void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);

FLINT_DLL void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A,
    const fmpq_mat_t B);

FLINT_DLL void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);

FLINT_DLL void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A,
    const fmpz_mat_t B);

FLINT_DLL void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A,
    const fmpq_mat_t B);

/* Trace *********************************************************************/

FLINT_DLL void fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat);

/* Determinant ***************************************************************/

FLINT_DLL void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat);

/* Nonsingular solving *******************************************************/

FLINT_DLL int fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A,
    const fmpq_mat_t B);

FLINT_DLL int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B);

int fmpq_mat_solve_fmpz_mat(fmpq_mat_t X, const fmpz_mat_t A,
    const fmpz_mat_t B);

/* Inverse *******************************************************************/

FLINT_DLL int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A);

/* Echelon form **************************************************************/

FLINT_DLL int fmpq_mat_pivot(slong * perm, fmpq_mat_t mat, slong r, slong c);

FLINT_DLL slong fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A);

FLINT_DLL slong fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A);

FLINT_DLL slong fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A);

/* Gram-Schmidt Orthogonalisation  *******************************************/

FLINT_DLL void fmpq_mat_gso(fmpq_mat_t B, const fmpq_mat_t A);

#ifdef __cplusplus
}
#endif

#endif

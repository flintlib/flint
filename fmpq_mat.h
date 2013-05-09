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
    len_t r;
    len_t c;
    fmpq ** rows;
} fmpq_mat_struct;

typedef fmpq_mat_struct fmpq_mat_t[1];

#define fmpq_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))
#define fmpq_mat_entry_num(mat,i,j) ((fmpz *)(&((*fmpq_mat_entry(mat,i,j)).num)))
#define fmpq_mat_entry_den(mat,i,j) ((fmpz *)(&((*fmpq_mat_entry(mat,i,j)).den)))

#define fmpq_mat_nrows(mat) ((mat)->r)
#define fmpq_mat_ncols(mat) ((mat)->c)

void fmpq_mat_init(fmpq_mat_t mat, len_t rows, len_t cols);

void fmpq_mat_clear(fmpq_mat_t mat);


void fmpq_mat_print(const fmpq_mat_t mat);

/* Random matrix generation **************************************************/

void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

/* Special matrices **********************************************************/

void fmpq_mat_hilbert_matrix(fmpq_mat_t mat);

/* Basic assignment **********************************************************/

void fmpq_mat_set(fmpq_mat_t dest, const fmpq_mat_t src);

void fmpq_mat_zero(fmpq_mat_t mat);

void fmpq_mat_one(fmpq_mat_t mat);

void fmpq_mat_transpose(fmpq_mat_t rop, const fmpq_mat_t op);

/* Addition, scalar multiplication  ******************************************/

void fmpq_mat_add(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2);

void fmpq_mat_sub(fmpq_mat_t mat, const fmpq_mat_t mat1, const fmpq_mat_t mat2);

void fmpq_mat_neg(fmpq_mat_t rop, const fmpq_mat_t op);

void fmpq_mat_scalar_mul_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x);

void fmpq_mat_scalar_div_fmpz(fmpq_mat_t rop, const fmpq_mat_t op, const fmpz_t x);

/* Basic comparison and properties *******************************************/

int fmpq_mat_equal(const fmpq_mat_t mat1, const fmpq_mat_t mat2);

int fmpq_mat_is_integral(const fmpq_mat_t mat);

int fmpq_mat_is_zero(const fmpq_mat_t mat);

static __inline__ int
fmpq_mat_is_empty(const fmpq_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
fmpq_mat_is_square(const fmpq_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Integer matrix conversion *************************************************/

int fmpq_mat_get_fmpz_mat(fmpz_mat_t dest, const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_entrywise(fmpz_mat_t num, fmpz_mat_t den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_matwise(fmpz_mat_t num, fmpz_t den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2,
        fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2);

void fmpq_mat_get_fmpz_mat_mod_fmpz(fmpz_mat_t dest, const fmpq_mat_t mat,
    const fmpz_t mod);

void fmpq_mat_set_fmpz_mat(fmpq_mat_t dest, const fmpz_mat_t src);

void fmpq_mat_set_fmpz_mat_div_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod,
    const fmpz_t div);

int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod,
    const fmpz_t mod);

/* Matrix multiplication *****************************************************/

void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);

void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A,
    const fmpq_mat_t B);

void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);

void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A,
    const fmpz_mat_t B);

void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A,
    const fmpq_mat_t B);

/* Trace *********************************************************************/

void fmpq_mat_trace(fmpq_t trace, const fmpq_mat_t mat);

/* Determinant ***************************************************************/

void fmpq_mat_det(fmpq_t det, const fmpq_mat_t mat);

/* Nonsingular solving *******************************************************/

int fmpq_mat_solve_fraction_free(fmpq_mat_t X, const fmpq_mat_t A,
    const fmpq_mat_t B);

int fmpq_mat_solve_dixon(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B);

/* Inverse *******************************************************************/

int fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A);

/* Echelon form **************************************************************/

int fmpq_mat_pivot(len_t * perm, fmpq_mat_t mat, len_t r, len_t c);

len_t fmpq_mat_rref_classical(fmpq_mat_t B, const fmpq_mat_t A);

len_t fmpq_mat_rref_fraction_free(fmpq_mat_t B, const fmpq_mat_t A);

len_t fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A);

#ifdef __cplusplus
}
#endif

#endif

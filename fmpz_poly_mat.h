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

#ifndef POLY_MAT_H
#define POLY_MAT_H

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "mat_common.h"

typedef struct
{
    fmpz_poly_struct * entries;
    long r;
    long c;
    fmpz_poly_struct ** rows;
} fmpz_poly_mat_struct;

typedef fmpz_poly_mat_struct fmpz_poly_mat_t[1];

#define fmpz_poly_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))

/* Memory management */

void fmpz_poly_mat_init(fmpz_poly_mat_t mat, long rows, long cols);

void fmpz_poly_mat_init_set(fmpz_poly_mat_t mat, const fmpz_poly_mat_t src);

void fmpz_poly_mat_swap(fmpz_poly_mat_t mat1, fmpz_poly_mat_t mat2);

void fmpz_poly_mat_set(fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2);

void fmpz_poly_mat_clear(fmpz_poly_mat_t mat);

/* Comparison */

int fmpz_poly_mat_equal(const fmpz_poly_mat_t mat1, const fmpz_poly_mat_t mat2);

int fmpz_poly_mat_is_zero(const fmpz_poly_mat_t mat);

static __inline__ int
fmpz_poly_mat_is_empty(const fmpz_poly_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
fmpz_poly_mat_is_square(const fmpz_poly_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Standard matrices */

void fmpz_poly_mat_zero(fmpz_poly_mat_t mat);

void fmpz_poly_mat_unit(fmpz_poly_mat_t mat);

/* Random matrices */

void fmpz_poly_mat_randtest(fmpz_poly_mat_t mat, flint_rand_t state, long len, mp_bitcnt_t bits);

/* Input and output */

void fmpz_poly_mat_print(const fmpz_poly_mat_t mat, const char * x);

/* Norms */

long fmpz_poly_mat_max_bits(const fmpz_poly_mat_t A);

long fmpz_poly_mat_max_length(const fmpz_poly_mat_t A);

/* Matrix arithmetic */

void fmpz_poly_mat_add(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_sub(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_neg(fmpz_poly_mat_t B, const fmpz_poly_mat_t A);

void fmpz_poly_mat_mul(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_mul_classical(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

void fmpz_poly_mat_mul_KS(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

/* Evaluation */

void fmpz_poly_mat_evaluate_fmpz(fmpz_mat_t B, const fmpz_poly_mat_t A, const fmpz_t x);

/* Row reduction */

int fmpz_poly_mat_pivot(long * perm, fmpz_poly_mat_t A, long r, long c);

long fmpz_poly_mat_rowreduce(long * perm, fmpz_poly_mat_t B, fmpz_poly_t den, const fmpz_poly_mat_t A, int options);

void fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A);

void fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A);

/* TBA

int fmpz_poly_mat_solve(fmpz_poly_mat_t X, fmpz_poly_t den, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

int fmpz_poly_mat_inv(fmpz_poly_mat_t B, fmpz_poly_t den, const fmpz_poly_mat_t A);

long fmpz_poly_mat_rank(const fmpz_poly_mat_t A);
*/

#endif

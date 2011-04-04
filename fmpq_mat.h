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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq.h"

typedef struct
{
    fmpq * entries;
    long r;
    long c;
    fmpq ** rows;
} fmpq_mat_struct;

typedef fmpq_mat_struct fmpq_mat_t[1];

#define fmpq_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))
#define fmpq_mat_entry_num(mat,i,j) ((fmpz *)(&((*fmpq_mat_entry(mat,i,j)).num)))
#define fmpq_mat_entry_den(mat,i,j) ((fmpz *)(&((*fmpq_mat_entry(mat,i,j)).den)))

void fmpq_mat_init(fmpq_mat_t mat, long rows, long cols);

void fmpq_mat_clear(fmpq_mat_t mat);

int fmpq_mat_equal(fmpq_mat_t mat1, fmpq_mat_t mat2);

void fmpq_mat_print(fmpq_mat_t mat);

void fmpq_mat_randbits(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

void fmpq_mat_randtest(fmpq_mat_t mat, flint_rand_t state, mp_bitcnt_t bits);

void fmpq_mat_get_fmpz_mat_rowwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_colwise(fmpz_mat_t num, fmpz * den,
    const fmpq_mat_t mat);

void fmpq_mat_get_fmpz_mat_rowwise_2(fmpz_mat_t num, fmpz_mat_t num2,
        fmpz * den, const fmpq_mat_t mat, const fmpq_mat_t mat2);


void fmpq_mat_mul_direct(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);
void fmpq_mat_mul_cleared(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);
void fmpq_mat_mul(fmpq_mat_t C, const fmpq_mat_t A, const fmpq_mat_t B);
void fmpq_mat_mul_fmpz_mat(fmpq_mat_t C, const fmpq_mat_t A, const fmpz_mat_t B);
void fmpq_mat_mul_r_fmpz_mat(fmpq_mat_t C, const fmpz_mat_t A, const fmpq_mat_t B);

void fmpq_mat_det(fmpq_t det, fmpq_mat_t mat);
void fmpq_mat_solve_mat(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B);

int fmpq_mat_set_fmpz_mat_mod_fmpz(fmpq_mat_t X, const fmpz_mat_t Xmod, const fmpz_t mod);

void fmpq_mat_inv(fmpq_mat_t B, const fmpq_mat_t A);
void fmpq_mat_rref(fmpq_mat_t B, const fmpq_mat_t A);

#endif

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
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#ifndef FMPZ_MAT_H
#define FMPZ_MAT_H

#include <mpir.h>
#include "fmpz.h"

typedef struct
{
    fmpz * entries;
    long r;
    long c;
    fmpz ** rows;
} fmpz_mat_struct;

/* fmpz_mat_t allows reference-like semantics for fmpz_mat_struct */
typedef fmpz_mat_struct fmpz_mat_t[1];

void fmpz_mat_init(fmpz_mat_t mat, long rows, long cols);

void fmpz_mat_clear(fmpz_mat_t mat);

void fmpz_mat_print(fmpz_mat_t mat); 

void fmpz_mat_print_pretty(fmpz_mat_t mat);

void fmpz_mat_randinit(fmpz_randstate_t state);

void fmpz_mat_randclear(fmpz_randstate_t state);

void fmpz_mat_randbits(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_mat_randtest(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_mat_randintrel(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits);

void fmpz_mat_randsimdioph(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits, mp_bitcnt_t bits2);

void fmpz_mat_randntrulike(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits, ulong q);

void fmpz_mat_randntrulike2(fmpz_mat_t mat, fmpz_randstate_t state, mp_bitcnt_t bits, ulong q);

void fmpz_mat_randajtai(fmpz_mat_t mat, fmpz_randstate_t state, double alpha);

void fmpz_mat_swap(fmpz_mat_t mat1, fmpz_mat_t mat2);


/* Linear algebra operations */

void fmpz_mat_mul(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

#define ROWREDUCE_DETONLY 0
#define ROWREDUCE_ECHELON_FORM 1
#define ROWREDUCE_REDUCED_ECHELON_FORM 2

long _fmpz_mat_rowreduce(fmpz ** a, long n, long m, int variant);

void _fmpz_mat_det_2x2(fmpz_t det, fmpz ** const x);
void _fmpz_mat_det_3x3(fmpz_t det, fmpz ** const x);
void _fmpz_mat_det_4x4(fmpz_t det, fmpz ** const x);
void _fmpz_mat_det_rowreduce(fmpz_t det, const fmpz_mat_t A);
void fmpz_mat_det(fmpz_t det, const fmpz_mat_t A);

long fmpz_mat_rank(const fmpz_mat_t A);

void _fmpz_mat_solve_2x2(fmpz * x, fmpz_t d, fmpz ** const a, const fmpz * b);
void _fmpz_mat_solve_3x3(fmpz * x, fmpz_t d, fmpz ** const a, const fmpz * b);
void _fmpz_mat_solve_rowreduce(fmpz * xnum, fmpz_t xden, const fmpz_mat_t A, const fmpz * b);
void fmpz_mat_solve(fmpz * xnum, fmpz_t xden, const fmpz_mat_t A, const fmpz * b);

void fmpz_mat_randrank(fmpz_mat_t mat, fmpz_randstate_t state, long rank, mp_bitcnt_t bits);
void fmpz_mat_randdet(fmpz_mat_t mat, fmpz_randstate_t state, const fmpz_t det);
void fmpz_mat_randops(fmpz_mat_t mat, fmpz_randstate_t state, long count);
int fmpz_mat_randpermdiag(fmpz_mat_t mat, fmpz_randstate_t state, const fmpz * diag, long n);


#endif


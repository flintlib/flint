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

#ifndef NMOD_MAT_H
#define NMOD_MAT_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#define ulong unsigned long

#include <mpir.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "mat_common.h"

typedef struct
{
    mp_limb_t * entries;
    long r;
    long c;
    mp_limb_t ** rows;
    nmod_t mod;
}
nmod_mat_struct;

/* fmpz_mat_t allows reference-like semantics for fmpz_mat_struct */
typedef nmod_mat_struct nmod_mat_t[1];

#define nmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])

static __inline__
void
_nmod_mat_set_mod(nmod_mat_t mat, mp_limb_t n)
{
    mat->mod.n = n;
    mat->mod.ninv = n_preinvert_limb(n);
    count_leading_zeros(mat->mod.norm, n);
}

/* Memory management */
void nmod_mat_init(nmod_mat_t mat, long rows, long cols, mp_limb_t n);
void nmod_mat_init_set(nmod_mat_t mat, const nmod_mat_t src);
void nmod_mat_clear(nmod_mat_t mat);

void nmod_mat_window_init(nmod_mat_t window, const nmod_mat_t mat, long r1, long c1, long r2, long c2);
void nmod_mat_window_clear(nmod_mat_t window);

/* Random matrix generation */
void nmod_mat_randtest(nmod_mat_t mat, flint_rand_t state);
void nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state);
int nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state, 
                 const mp_limb_t * diag, long n);
void nmod_mat_randrank(nmod_mat_t, flint_rand_t state, long rank);
void nmod_mat_randops(nmod_mat_t mat, long count, flint_rand_t state);
void nmod_mat_randtril(nmod_mat_t mat, flint_rand_t state, int unit);


void nmod_mat_print_pretty(const nmod_mat_t mat);

int nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2);


void nmod_mat_set(nmod_mat_t B, const nmod_mat_t A);
void nmod_mat_transpose(nmod_mat_t B, const nmod_mat_t A);

/* Arithmetic */
void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_mul_strassen(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);


int _nmod_mat_pivot(mp_limb_t ** rows, long n, long start_row, long col);

long _nmod_mat_rowreduce_1(nmod_mat_t mat, int options);
long _nmod_mat_rowreduce_2(nmod_mat_t mat, int options);
long _nmod_mat_rowreduce_r(nmod_mat_t mat, int options);
long _nmod_mat_rowreduce(nmod_mat_t mat, int options);

mp_limb_t _nmod_mat_fast_rowreduce_modulus_1(long rows, long cols, int proved);
mp_limb_t _nmod_mat_fast_rowreduce_modulus_2(long rows, long cols, int proved);
mp_limb_t _nmod_mat_fast_rowreduce_modulus(long rows, long cols, int proved);


mp_limb_t _nmod_mat_det_rowreduce(nmod_mat_t A);
mp_limb_t nmod_mat_det(const nmod_mat_t A);

long nmod_mat_rank(const nmod_mat_t A);

void _nmod_mat_solve_lu_precomp(mp_limb_t * b, mp_limb_t ** const LU, long n, nmod_t mod);
int nmod_mat_solve(mp_limb_t * x, const nmod_mat_t A, const mp_limb_t * b);
int nmod_mat_solve_mat(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B);

int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A);

/* New solvers */

void nmod_mat_solve_tril(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
void nmod_mat_solve_tril_recursive(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);
void nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L, const nmod_mat_t B, int unit);

#endif

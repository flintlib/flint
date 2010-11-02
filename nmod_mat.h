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

#include <stdlib.h>
#include <mpir.h>
#include "longlong.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

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

/* Memory management */
void nmod_mat_init(nmod_mat_t mat, long rows, long cols, mp_limb_t n);
void nmod_mat_clear(nmod_mat_t mat);

void nmod_mat_randtest(nmod_mat_t mat);

void nmod_mat_print_pretty(nmod_mat_t mat);

int nmod_mat_equal(const nmod_mat_t mat1, const nmod_mat_t mat2);

/* Arithmetic */
void nmod_mat_add(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void nmod_mat_mul(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

void _nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void _nmod_mat_mul_classical_1(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void _nmod_mat_mul_classical_1r(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long run_length);
void _nmod_mat_mul_classical_2(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);
void _nmod_mat_mul_classical_2r(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long run_length);
void _nmod_mat_mul_classical_3(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

void _nmod_mat_mul_blocked_1(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);
void _nmod_mat_mul_blocked_1b(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);
void _nmod_mat_mul_blocked_2(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);
void _nmod_mat_mul_blocked_2b(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);
void _nmod_mat_mul_blocked_3(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);
void _nmod_mat_mul_blocked_3b(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size);

#endif

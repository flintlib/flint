/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef BOOL_MAT_H
#define BOOL_MAT_H

#ifdef BOOL_MAT_INLINES_C
#define BOOL_MAT_INLINE
#else
#define BOOL_MAT_INLINE static inline
#endif

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int *entries;
    slong r;
    slong c;
    int **rows;
}
bool_mat_struct;

typedef bool_mat_struct bool_mat_t[1];

#define bool_mat_nrows(mat) ((mat)->r)
#define bool_mat_ncols(mat) ((mat)->c)

BOOL_MAT_INLINE int
bool_mat_get_entry(const bool_mat_t mat, slong i, slong j)
{
    return mat->rows[i][j];
}

BOOL_MAT_INLINE void
bool_mat_set_entry(bool_mat_t mat, slong i, slong j, int value)
{
    mat->rows[i][j] = value;
}

/* Memory management */

void bool_mat_init(bool_mat_t mat, slong r, slong c);

void bool_mat_clear(bool_mat_t mat);

BOOL_MAT_INLINE void
bool_mat_swap(bool_mat_t mat1, bool_mat_t mat2)
{
    bool_mat_struct t = *mat1;
    *mat1 = *mat2;
    *mat2 = t;
}

/* Conversions */

void bool_mat_set(bool_mat_t dest, const bool_mat_t src);

/* Random generation */

void bool_mat_randtest(bool_mat_t mat, flint_rand_t state);

void bool_mat_randtest_diagonal(bool_mat_t mat, flint_rand_t state);

void bool_mat_randtest_nilpotent(bool_mat_t mat, flint_rand_t state);

/* I/O */

#ifdef FLINT_HAVE_FILE
void bool_mat_fprint(FILE * file, const bool_mat_t mat);
#endif

void bool_mat_print(const bool_mat_t mat);

/* Comparisons */

int bool_mat_equal(const bool_mat_t mat1, const bool_mat_t mat2);

int bool_mat_any(const bool_mat_t mat);

int bool_mat_all(const bool_mat_t mat);

int bool_mat_is_diagonal(const bool_mat_t mat);

int bool_mat_is_lower_triangular(const bool_mat_t mat);

int bool_mat_is_transitive(const bool_mat_t mat);

int bool_mat_is_nilpotent(const bool_mat_t mat);

BOOL_MAT_INLINE int
bool_mat_is_empty(const bool_mat_t mat)
{
    return (mat->r == 0) || (mat->c == 0);
}

BOOL_MAT_INLINE int
bool_mat_is_square(const bool_mat_t mat)
{
    return (mat->r == mat->c);
}

/* Special matrices */

void bool_mat_zero(bool_mat_t mat);

void bool_mat_one(bool_mat_t mat);

void bool_mat_directed_path(bool_mat_t mat);

void bool_mat_directed_cycle(bool_mat_t mat);

/* Transpose */

void bool_mat_transpose(bool_mat_t mat1, const bool_mat_t mat2);

/* Arithmetic */

void bool_mat_complement(bool_mat_t mat1, const bool_mat_t mat2);

void bool_mat_add(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2);

void bool_mat_mul(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2);

void bool_mat_mul_entrywise(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2);

void bool_mat_pow_ui(bool_mat_t B, const bool_mat_t A, ulong exp);

BOOL_MAT_INLINE void
bool_mat_sqr(bool_mat_t B, const bool_mat_t A)
{
    bool_mat_mul(B, A, A);
}

/* Special functions */

int bool_mat_trace(const bool_mat_t mat);

slong bool_mat_nilpotency_degree(const bool_mat_t mat);

void bool_mat_transitive_closure(bool_mat_t dest, const bool_mat_t src);

slong bool_mat_get_strongly_connected_components(slong *partition, const bool_mat_t A);

slong bool_mat_all_pairs_longest_walk(fmpz_mat_t B, const bool_mat_t A);

#ifdef __cplusplus
}
#endif

#endif

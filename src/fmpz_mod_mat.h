/*
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MAT_H
#define FMPZ_MOD_MAT_H

#ifdef FMPZ_MOD_MAT_INLINES_C
#define FMPZ_MOD_MAT_INLINE
#else
#define FMPZ_MOD_MAT_INLINE static inline
#endif

#include "fmpz_mat.h"
#include "fmpz_mod_types.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF 10

/* Element access  ********************************************************/

FMPZ_MOD_MAT_INLINE
fmpz * fmpz_mod_mat_entry(const fmpz_mod_mat_t mat, slong i, slong j)
{
    return fmpz_mat_entry(mat, i, j);
}

void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val, const fmpz_mod_ctx_t ctx);
void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j, const fmpz_mod_ctx_t ctx);

/* Memory management  ********************************************************/

void fmpz_mod_mat_init(fmpz_mod_mat_t mat, slong rows, slong cols, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_init_set(fmpz_mod_mat_t mat, const fmpz_mod_mat_t src, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_clear(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

/* Basic manipulation  ********************************************************/

FMPZ_MOD_MAT_INLINE
slong fmpz_mod_mat_nrows(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_nrows(mat);
}

FMPZ_MOD_MAT_INLINE
slong fmpz_mod_mat_ncols(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_ncols(mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_one(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    if (*ctx->n == WORD(1))
        fmpz_mat_zero(mat);
    else
        fmpz_mat_one(mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_zero(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_zero(mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_empty(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_is_empty(mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_square(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    return fmpz_mat_is_square(mat);
}

void fmpz_mod_mat_swap(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx);
void fmpz_mod_mat_swap_entrywise(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

/* Conversions */

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_set_nmod_mat(fmpz_mod_mat_t A, const nmod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_set_nmod_mat_unsigned(A, B);
}

void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t A, const fmpz_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MAT_INLINE void _fmpz_mod_mat_reduce(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_mat_set_fmpz_mat(mat, mat, ctx);
}

/* Random matrix generation */

void fmpz_mod_mat_randtest(fmpz_mod_mat_t mat, flint_rand_t state, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state, slong rank, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_randtril(fmpz_mod_mat_t mat, flint_rand_t state, int unit, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_randtriu(fmpz_mod_mat_t mat, flint_rand_t state, int unit, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_randops(fmpz_mod_mat_t mat, flint_rand_t state, slong count, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_randops(mat, state, count);
    _fmpz_mod_mat_reduce(mat, ctx);
}

/* Windows and concatenation */

void fmpz_mod_mat_window_init(fmpz_mod_mat_t window, const fmpz_mod_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_window_clear(fmpz_mod_mat_t window, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_concat_horizontal(fmpz_mod_mat_t res,
                                    const fmpz_mod_mat_t mat1,
                                    const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_concat_horizontal(res, mat1, mat2);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_concat_vertical(fmpz_mod_mat_t res,
                                  const fmpz_mod_mat_t mat1,
                                  const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_concat_vertical(res, mat1, mat2);
}

/* Input/output */

#ifdef FLINT_HAVE_FILE
int fmpz_mod_mat_fprint(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);
int fmpz_mod_mat_fprint_pretty(FILE * file, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);
#endif

int fmpz_mod_mat_print(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);
void fmpz_mod_mat_print_pretty(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

/* Comparison */

int fmpz_mod_mat_equal(const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2, const fmpz_mod_ctx_t ctx);
int fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);
int fmpz_mod_mat_is_one(const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

/* Transpose */

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_transpose(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_transpose(B, A);
}

/* Addition and subtraction */

void fmpz_mod_mat_add(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_sub(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_neg(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

/* Matrix-scalar arithmetic */

void fmpz_mod_mat_scalar_mul_si(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, slong c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_scalar_mul_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong c, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_scalar_mul_fmpz(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, fmpz_t c, const fmpz_mod_ctx_t ctx);

/* Matrix multiplication */

void fmpz_mod_mat_mul(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_mat_mul_classical_threaded_pool_op(fmpz_mod_mat_t D,
      const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B,
                      int op, thread_pool_handle * threads, slong num_threads, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C,
                               const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_sqr(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_submul(fmpz_mod_mat_t D, const fmpz_mod_mat_t C,
                               const fmpz_mod_mat_t A, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_mul_fmpz_vec(fmpz * c, const fmpz_mod_mat_t A,
                                                   const fmpz * b, slong blen, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_mul_fmpz_vec_ptr(fmpz * const * c,
                   const fmpz_mod_mat_t A, const fmpz * const * b, slong blen, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen,
                                                       const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_fmpz_vec_mul_ptr(fmpz * const * c,
                   const fmpz * const * a, slong alen, const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

/* Trace */

void fmpz_mod_mat_trace(fmpz_t trace, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

/* Gaussian elimination *********************************************/

void fmpz_mod_mat_det(fmpz_t res, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

slong fmpz_mod_mat_rref(fmpz_mod_mat_t res, const fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

slong fmpz_mod_mat_reduce_row(fmpz_mod_mat_t A, slong * P, slong * L, slong m, const fmpz_mod_ctx_t ctx);

slong fmpz_mod_mat_lu(slong * P, fmpz_mod_mat_t A, int rank_check, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_solve_triu(fmpz_mod_mat_t X, const fmpz_mod_mat_t L,
                                             const fmpz_mod_mat_t B, int unit, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_solve_tril(fmpz_mod_mat_t X, const fmpz_mod_mat_t L,
                                             const fmpz_mod_mat_t B, int unit, const fmpz_mod_ctx_t ctx);

int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                                                       const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

int fmpz_mod_mat_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                                                       const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx);

int fmpz_mod_mat_inv(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

slong fmpz_mod_mat_rank(const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

slong fmpz_mod_mat_nullspace(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx);

/* Howell and strong echelon form ***********************************/

slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

void fmpz_mod_mat_strong_echelon_form(fmpz_mod_mat_t mat, const fmpz_mod_ctx_t ctx);

/* Transforms ****************************************************************/

void fmpz_mod_mat_similarity(fmpz_mod_mat_t A, slong r, fmpz_t d, const fmpz_mod_ctx_t ctx);

/* Permutations ************************************************************/

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_swap_rows(fmpz_mod_mat_t mat, slong * perm, slong r, slong s, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_swap_rows(mat, perm, r, s);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_invert_rows(fmpz_mod_mat_t mat, slong * perm, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_invert_rows(mat, perm);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_swap_cols(fmpz_mod_mat_t mat, slong * perm, slong r, slong s, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_swap_cols(mat, perm, r, s);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_invert_cols(fmpz_mod_mat_t mat, slong * perm, const fmpz_mod_ctx_t ctx)
{
    fmpz_mat_invert_cols(mat, perm);
}

#ifdef __cplusplus
}
#endif

#endif

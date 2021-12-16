/*
    Copyright (C) 2017 Luca De Feo
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_MAT_H
#define FMPZ_MOD_MAT_H

#ifdef FMPZ_MOD_MAT_INLINES_C
#define FMPZ_MOD_MAT_INLINE FLINT_DLL
#else
#define FMPZ_MOD_MAT_INLINE static __inline__
#endif

#include "flint.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF 10
#define FMPZ_MOD_MAT_LU_RECURSIVE_CUTOFF 4
#define FMPZ_MOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define FMPZ_MOD_MAT_SOLVE_TRI_COLS_CUTOFF 64

typedef struct
{
    fmpz_mat_t mat;
    fmpz_t mod;
}
fmpz_mod_mat_struct;

/* fmpz_mod_mat_t allows reference-like semantics for fmpz_mod_mat_struct */
typedef fmpz_mod_mat_struct fmpz_mod_mat_t[1];

/* Element access  ********************************************************/

FMPZ_MOD_MAT_INLINE
fmpz * fmpz_mod_mat_entry(const fmpz_mod_mat_t mat, slong i, slong j)
{
    return fmpz_mat_entry(mat->mat, i, j);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_set_entry(fmpz_mod_mat_t mat, slong i, slong j, const fmpz_t val)
{
    fmpz_set(fmpz_mat_entry(mat->mat, i, j), val);
}

/* Memory management  ********************************************************/

FLINT_DLL void fmpz_mod_mat_init(fmpz_mod_mat_t mat, slong rows, slong cols, const fmpz_t n);

FLINT_DLL void fmpz_mod_mat_init_set(fmpz_mod_mat_t mat, const fmpz_mod_mat_t src);

FLINT_DLL void fmpz_mod_mat_clear(fmpz_mod_mat_t mat);

/* Basic manipulation  ********************************************************/

FMPZ_MOD_MAT_INLINE
slong fmpz_mod_mat_nrows(const fmpz_mod_mat_t mat)
{
    return fmpz_mat_nrows(mat->mat);
}

FMPZ_MOD_MAT_INLINE
slong fmpz_mod_mat_ncols(const fmpz_mod_mat_t mat)
{
    return fmpz_mat_ncols(mat->mat);
}

FMPZ_MOD_MAT_INLINE
void _fmpz_mod_mat_set_mod(fmpz_mod_mat_t mat, const fmpz_t n)
{
    fmpz_set(mat->mod, n);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_one(fmpz_mod_mat_t mat)
{
    fmpz_mat_one(mat->mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_zero(fmpz_mod_mat_t mat)
{
    fmpz_mat_zero(mat->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_empty(const fmpz_mod_mat_t mat)
{
	    return fmpz_mat_is_empty(mat->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_square(const fmpz_mod_mat_t mat)
{
	    return fmpz_mat_is_square(mat->mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_swap(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
{
    if (mat1 != mat2)
    {
        fmpz_mod_mat_struct tmp;

        tmp = *mat1;
        *mat1 = *mat2;
        *mat2 = tmp;
    }
}

FMPZ_MOD_MAT_INLINE void
fmpz_mod_mat_swap_entrywise(fmpz_mod_mat_t mat1, fmpz_mod_mat_t mat2)
{
    slong i, j;

    for (i = 0; i < fmpz_mod_mat_nrows(mat1); i++)
        for (j = 0; j < fmpz_mod_mat_ncols(mat1); j++)
            fmpz_swap(fmpz_mod_mat_entry(mat2, i, j), fmpz_mod_mat_entry(mat1, i, j));
}

FLINT_DLL void _fmpz_mod_mat_reduce(fmpz_mod_mat_t mat);

/* Random matrix generation */
FLINT_DLL void fmpz_mod_mat_randtest(fmpz_mod_mat_t mat, flint_rand_t state);

FLINT_DLL void fmpz_mod_mat_randrank(fmpz_mod_mat_t mat, flint_rand_t state,
                                                                   slong rank);

FLINT_DLL void fmpz_mod_mat_randtril(fmpz_mod_mat_t mat, flint_rand_t state,
                                                                     int unit);

FLINT_DLL void fmpz_mod_mat_randtriu(fmpz_mod_mat_t mat, flint_rand_t state,
                                                                     int unit);

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_randops(fmpz_mod_mat_t mat, slong count, flint_rand_t state)
{
    fmpz_mat_randops(mat->mat, state, count);
    _fmpz_mod_mat_reduce(mat);
}

/* Windows and concatenation */

FLINT_DLL void fmpz_mod_mat_window_init(fmpz_mod_mat_t window, const fmpz_mod_mat_t mat,
                              slong r1, slong c1, slong r2, slong c2);

FLINT_DLL void fmpz_mod_mat_window_clear(fmpz_mod_mat_t window);

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_concat_horizontal(fmpz_mod_mat_t res,
                                    const fmpz_mod_mat_t mat1,
                                    const fmpz_mod_mat_t mat2)
{
    fmpz_mat_concat_horizontal(res->mat, mat1->mat, mat2->mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_concat_vertical(fmpz_mod_mat_t res,
                                  const fmpz_mod_mat_t mat1,
                                  const fmpz_mod_mat_t mat2)
{
    fmpz_mat_concat_vertical(res->mat, mat1->mat, mat2->mat);
}

/* Input/output */
FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_fprint(FILE * file, const fmpz_mod_mat_t mat)
{
    return fmpz_mat_fprint(file, mat->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_fprint_pretty(FILE * file, const fmpz_mod_mat_t mat)
{
    return fmpz_mat_fprint_pretty(file, mat->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_print(const fmpz_mod_mat_t mat)
{
    return fmpz_mat_print(mat->mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_print_pretty(const fmpz_mod_mat_t mat)
{
    fmpz_mat_print_pretty(mat->mat);
}

/* Comparison */
FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_equal(const fmpz_mod_mat_t mat1, const fmpz_mod_mat_t mat2)
{
    return fmpz_equal(mat1->mod, mat2->mod) && fmpz_mat_equal(mat1->mat, mat2->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_zero(const fmpz_mod_mat_t mat)
{
    return fmpz_mat_is_zero(mat->mat);
}

FMPZ_MOD_MAT_INLINE
int fmpz_mod_mat_is_one(const fmpz_mod_mat_t mat)
{
    return fmpz_is_one(mat->mod) || fmpz_mat_is_one(mat->mat);
}


/* Set and transpose */
FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_set(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
{
    fmpz_set(B->mod, A->mod);
    fmpz_mat_set(B->mat, A->mat);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_transpose(fmpz_mod_mat_t B, const fmpz_mod_mat_t A)
{
    fmpz_mat_transpose(B->mat, A->mat);
}

/* Conversions */

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_set_nmod_mat(fmpz_mod_mat_t A, const nmod_mat_t B)
{
    fmpz_mat_set_nmod_mat_unsigned(A->mat, B);
}

FLINT_DLL void fmpz_mod_mat_set_fmpz_mat(fmpz_mod_mat_t A, const fmpz_mat_t B);

FLINT_DLL void fmpz_mod_mat_get_fmpz_mat(fmpz_mat_t A, const fmpz_mod_mat_t B);

/* Addition and subtraction */

FLINT_DLL void fmpz_mod_mat_add(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B);

FLINT_DLL void fmpz_mod_mat_sub(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B);

FLINT_DLL void fmpz_mod_mat_neg(fmpz_mod_mat_t B, const fmpz_mod_mat_t A);

/* Matrix-scalar arithmetic */

FLINT_DLL void fmpz_mod_mat_scalar_mul_si(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, slong c);

FLINT_DLL void fmpz_mod_mat_scalar_mul_ui(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, ulong c);

FLINT_DLL void fmpz_mod_mat_scalar_mul_fmpz(fmpz_mod_mat_t B, const fmpz_mod_mat_t A, fmpz_t c);

/* Matrix multiplication */

FLINT_DLL void fmpz_mod_mat_mul(fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B);

FLINT_DLL void _fmpz_mod_mat_mul_classical_threaded_pool_op(fmpz_mod_mat_t D,
      const fmpz_mod_mat_t C, const fmpz_mod_mat_t A, const fmpz_mod_mat_t B,
                      int op, thread_pool_handle * threads, slong num_threads);

FLINT_DLL void fmpz_mod_mat_mul_classical_threaded(fmpz_mod_mat_t C,
                               const fmpz_mod_mat_t A, const fmpz_mod_mat_t B);

FLINT_DLL void fmpz_mod_mat_sqr(fmpz_mod_mat_t B, const fmpz_mod_mat_t A);

FLINT_DLL void fmpz_mod_mat_submul(fmpz_mod_mat_t D, const fmpz_mod_mat_t C,
                               const fmpz_mod_mat_t A, const fmpz_mod_mat_t B);

FLINT_DLL void fmpz_mod_mat_mul_fmpz_vec(fmpz * c, const fmpz_mod_mat_t A,
                                                   const fmpz * b, slong blen);

FLINT_DLL void fmpz_mod_mat_mul_fmpz_vec_ptr(fmpz * const * c,
                   const fmpz_mod_mat_t A, const fmpz * const * b, slong blen);

FLINT_DLL void fmpz_mod_mat_fmpz_vec_mul(fmpz * c, const fmpz * a, slong alen,
                                                       const fmpz_mod_mat_t B);

FLINT_DLL void fmpz_mod_mat_fmpz_vec_mul_ptr(fmpz * const * c,
                   const fmpz * const * a, slong alen, const fmpz_mod_mat_t B);

/* Trace */

FLINT_DLL void fmpz_mod_mat_trace(fmpz_t trace, const fmpz_mod_mat_t mat);

/* Gaussian elimination *********************************************/

FLINT_DLL slong fmpz_mod_mat_rref(slong * perm, fmpz_mod_mat_t mat);

FLINT_DLL slong _fmpz_mod_mat_reduce_row(fmpz_mod_mat_t A, slong * P, slong * L, 
                                            slong m, const fmpz_mod_ctx_t ctx);

FLINT_DLL slong fmpz_mod_mat_reduce_row(fmpz_mod_mat_t A, slong * P, slong * L,
                                                                      slong m);

FLINT_DLL slong fmpz_mod_mat_lu(slong * P, fmpz_mod_mat_t A, int rank_check);

FLINT_DLL slong fmpz_mod_mat_lu_classical(slong * P, fmpz_mod_mat_t A,
                                                               int rank_check);

FLINT_DLL slong fmpz_mod_mat_lu_recursive(slong * P, fmpz_mod_mat_t A,
                                                               int rank_check);

FLINT_DLL void fmpz_mod_mat_solve_triu(fmpz_mod_mat_t X, const fmpz_mod_mat_t L,
                                             const fmpz_mod_mat_t B, int unit);

FLINT_DLL void fmpz_mod_mat_solve_triu_classical(fmpz_mod_mat_t X,
                    const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit);

FLINT_DLL void fmpz_mod_mat_solve_triu_recursive(fmpz_mod_mat_t X,
                    const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit);

FLINT_DLL void fmpz_mod_mat_solve_tril(fmpz_mod_mat_t X, const fmpz_mod_mat_t L,
                                             const fmpz_mod_mat_t B, int unit);

FLINT_DLL void fmpz_mod_mat_solve_tril_classical(fmpz_mod_mat_t X,
                    const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit);

FLINT_DLL void fmpz_mod_mat_solve_tril_recursive(fmpz_mod_mat_t X,
                    const fmpz_mod_mat_t L, const fmpz_mod_mat_t B, int unit);

FLINT_DLL int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                                                       const fmpz_mod_mat_t B);

FLINT_DLL int fmpz_mod_mat_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                                                       const fmpz_mod_mat_t B);

FLINT_DLL int fmpz_mod_mat_inv(fmpz_mod_mat_t B, fmpz_mod_mat_t A);

FLINT_DLL slong fmpz_mod_mat_rank(const fmpz_mod_mat_t A);

FLINT_DLL slong fmpz_mod_mat_nullspace(fmpz_mod_mat_t X, const fmpz_mod_mat_t A);

/* Howell and strong echelon form ***********************************/

FLINT_DLL slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat);

FLINT_DLL void fmpz_mod_mat_strong_echelon_form(fmpz_mod_mat_t mat);

/* Transforms ****************************************************************/

FLINT_DLL void fmpz_mod_mat_similarity(fmpz_mod_mat_t A, slong r, fmpz_t d);

/* Permutations ************************************************************/

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_swap_rows(fmpz_mod_mat_t mat, slong * perm, slong r, slong s)
{
    fmpz_mat_swap_rows(mat->mat, perm, r, s);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_invert_rows(fmpz_mod_mat_t mat, slong * perm)
{
    fmpz_mat_invert_rows(mat->mat, perm);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_swap_cols(fmpz_mod_mat_t mat, slong * perm, slong r, slong s)
{
    fmpz_mat_swap_cols(mat->mat, perm, r, s);
}

FMPZ_MOD_MAT_INLINE
void fmpz_mod_mat_invert_cols(fmpz_mod_mat_t mat, slong * perm)
{
    fmpz_mat_invert_cols(mat->mat, perm);
}

/* Inlines *******************************************************************/

FLINT_DLL void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j);

#ifdef __cplusplus
}
#endif

#endif


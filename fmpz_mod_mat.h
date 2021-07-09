/*
    Copyright (C) 2017 Luca De Feo

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

#ifdef __cplusplus
 extern "C" {
#endif

#define FMPZ_MOD_MAT_MUL_TRANSPOSE_CUTOFF 10

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

/* Howell and strong echelon form ***********************************/

FLINT_DLL slong fmpz_mod_mat_howell_form(fmpz_mod_mat_t mat);

FLINT_DLL void fmpz_mod_mat_strong_echelon_form(fmpz_mod_mat_t mat);

/* Inlines *******************************************************************/

FLINT_DLL void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j);

#ifdef __cplusplus
}
#endif

#endif


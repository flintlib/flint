/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "nmod_types.h"
#include "fmpz_mod_types.h"
#include "templates.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Memory management  ********************************************************/

void TEMPLATE(T, mat_init)(TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_init_set)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, mat_t) src,
                          const TEMPLATE(T, ctx_t) ctx);

FQ_MAT_TEMPLATES_INLINE slong
TEMPLATE(T, mat_nrows)(const TEMPLATE(T, mat_t) mat ,
                       const TEMPLATE(T, ctx_t) ctx)
{
    return mat->r;
}

FQ_MAT_TEMPLATES_INLINE slong
TEMPLATE(T, mat_ncols)(const TEMPLATE(T, mat_t) mat,
                       const TEMPLATE(T, ctx_t) ctx)
{
    return mat->c;
}

FQ_MAT_TEMPLATES_INLINE TEMPLATE(T, struct) *
TEMPLATE(T, mat_entry)(const TEMPLATE(T, mat_t) mat, slong i, slong j)
{
    return mat->rows[i] + j;
}

void TEMPLATE(T, mat_entry_set)(TEMPLATE(T, mat_t) mat, slong i, slong j,
                           const TEMPLATE(T, t) x,
                           const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_swap)(TEMPLATE(T, mat_t) mat1, TEMPLATE(T, mat_t) mat2,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_swap_entrywise)(TEMPLATE(T, mat_t) mat1,
		         TEMPLATE(T, mat_t) mat2, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_set)(TEMPLATE(T, mat_t) mat1, const TEMPLATE(T, mat_t) mat2,
                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_clear)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, mat_equal)(const TEMPLATE(T, mat_t) mat1,
                       const TEMPLATE(T, mat_t) mat2,
                       const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, mat_is_zero)(const TEMPLATE(T, mat_t) mat,
                         const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, mat_is_one)(const TEMPLATE(T, mat_t) mat,
		                         const TEMPLATE(T, ctx_t) ctx);

FQ_MAT_TEMPLATES_INLINE int
TEMPLATE(T, mat_is_empty)(const TEMPLATE(T, mat_t) mat,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (mat->r == 0) || (mat->c == 0);
}

FQ_MAT_TEMPLATES_INLINE int
TEMPLATE(T, mat_is_square)(const TEMPLATE(T, mat_t) mat,
                           const TEMPLATE(T, ctx_t) ctx)
{
    return (mat->r == mat->c);
}

FQ_MAT_TEMPLATES_INLINE void
TEMPLATE(T, mat_swap_rows)(TEMPLATE(T, mat_t) mat, slong * perm, slong r, slong s, const TEMPLATE(T, ctx_t) ctx)
{
    if (r != s && !TEMPLATE(T, mat_is_empty)(mat, ctx))
    {
        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        FLINT_SWAP(TEMPLATE(T, struct) *, mat->rows[r], mat->rows[s]);
    }
}

FQ_MAT_TEMPLATES_INLINE void
TEMPLATE(T, mat_invert_rows)(TEMPLATE(T, mat_t) mat, slong * perm, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < mat->r/2; i++)
        TEMPLATE(T, mat_swap_rows)(mat, perm, i, mat->r - i - 1, ctx);
}

void TEMPLATE(T, mat_swap_cols)(TEMPLATE(T, mat_t) mat, slong * perm, slong r, slong s, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_invert_cols)(TEMPLATE(T, mat_t) mat, slong * perm, const TEMPLATE(T, ctx_t) ctx);

/* Assignment  ***************************************************************/

void TEMPLATE(T, mat_zero)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_one)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

/* Conversions ***************************************************************/

void TEMPLATE(T, mat_set_nmod_mat) (TEMPLATE(T, mat_t) mat1,
                          const nmod_mat_t mat2, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_set_fmpz_mod_mat) (TEMPLATE(T, mat_t) mat1,
                      const fmpz_mod_mat_t mat2, const TEMPLATE(T, ctx_t) ctx);

/* Windows and concatenation */

void TEMPLATE(T, mat_window_init)(TEMPLATE(T, mat_t) window,
                             const TEMPLATE(T, mat_t) mat,
                             slong r1, slong c1, slong r2, slong c2,
                             const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_window_clear)(TEMPLATE(T, mat_t) window,
                              const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_concat_horizontal)(TEMPLATE(T, mat_t) res,
                           const TEMPLATE(T, mat_t) mat1,  const TEMPLATE(T, mat_t) mat2,
                           const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_concat_vertical)(TEMPLATE(T, mat_t) res,
                           const TEMPLATE(T, mat_t) mat1,  const TEMPLATE(T, mat_t) mat2,
                           const TEMPLATE(T, ctx_t) ctx);


/* Input and output  *********************************************************/

#ifdef FLINT_HAVE_FILE
int TEMPLATE(T, mat_fprint)(FILE * file, const TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, mat_fprint_pretty)(FILE * file, const TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);
#endif

int TEMPLATE(T, mat_print)(const TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, mat_print_pretty)(const TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);

/* TODO: Read functions */

/* Random matrix generation  *************************************************/

void TEMPLATE(T, mat_randtest)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_randrank)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          slong rank, const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, mat_randpermdiag)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                              TEMPLATE(T, struct) * diag, slong n,
                              const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_randops)(TEMPLATE(T, mat_t) mat, slong count,
                         flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_randtril)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          int unit, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_randtriu)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          int unit, const TEMPLATE(T, ctx_t) ctx);

/* Norms */

/* Transpose */

/* Addition and subtraction */

void TEMPLATE(T, mat_add)(TEMPLATE(T, mat_t) C,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_sub)(TEMPLATE(T, mat_t) C,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_neg)(TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_submul)(TEMPLATE(T, mat_t) D,
                        const TEMPLATE(T, mat_t) C,
                        const TEMPLATE(T, mat_t) A,
                        const TEMPLATE(T, mat_t) B,
                        const TEMPLATE(T, ctx_t) ctx);

/* Scalar operations */

/* Multiplication */

void TEMPLATE(T, mat_mul)(TEMPLATE(T, mat_t) C,
                     const TEMPLATE(T, mat_t) A,
                     const TEMPLATE(T, mat_t) B,
                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_mul_classical)(TEMPLATE(T, mat_t) C,
                               const TEMPLATE(T, mat_t) A,
                               const TEMPLATE(T, mat_t) B,
                               const TEMPLATE(T, ctx_t) ctx);


void TEMPLATE(T, mat_mul_KS)(TEMPLATE(T, mat_t) C,
                        const TEMPLATE(T, mat_t) A,
                        const TEMPLATE(T, mat_t) B,
                        const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, mat_lu)(slong * P,
                    TEMPLATE(T, mat_t) A,
                    int rank_check,
                    const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, mat_lu_recursive)(slong * P,
                              TEMPLATE(T, mat_t) A,
                              int rank_check,
                              const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, mat_lu_classical)(slong * P, TEMPLATE(T, mat_t) A, int rank_check,
                              const TEMPLATE(T, ctx_t) ctx);

/* Inverse *******************************************************************/

int TEMPLATE(T, mat_inv)(TEMPLATE(T, mat_t) B, TEMPLATE(T, mat_t) A,
                                   const TEMPLATE(T, ctx_t) ctx);

/* Solving *******************************************************************/

slong TEMPLATE(T, mat_rref)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);


slong TEMPLATE(T, mat_reduce_row)(TEMPLATE(T, mat_t) A, slong * P, slong * L,
                                         slong m, const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, mat_nullspace)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) A,
                           const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, mat_rank)(const TEMPLATE(T, mat_t) A,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_solve_tril)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) L,
                            const TEMPLATE(T, mat_t) B, int unit,
                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_solve_tril_classical)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) L,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_solve_tril_recursive)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) L,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_solve_triu)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) U,
                            const TEMPLATE(T, mat_t) B, int unit,
                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_solve_triu_classical)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) U,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE(T, mat_solve_triu_recursive)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) U,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_mul_vec)(TEMPLATE(T, struct) * c,
                                    const TEMPLATE(T, mat_t) A,
                                    const TEMPLATE(T, struct) * b, slong blen,
                                    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_mul_vec_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, mat_t) A,
                            const TEMPLATE(T, struct) * const * b, slong blen,
                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_vec_mul)(TEMPLATE(T, struct) * c,
                                    const TEMPLATE(T, struct) * a, slong alen,
                                    const TEMPLATE(T, mat_t) B,
                                    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_vec_mul_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, struct) * const * a, slong alen,
                            const TEMPLATE(T, mat_t) B,
                            const TEMPLATE(T, ctx_t) ctx);

/* Nonsingular solving *******************************************************/

int TEMPLATE(T, mat_solve)(TEMPLATE(T, mat_t) X,
               const TEMPLATE(T, mat_t A), const TEMPLATE(T, mat_t) C,
                                                 const TEMPLATE(T, ctx_t) ctx);

/* Solving *******************************************************************/

int TEMPLATE(T, mat_can_solve)(TEMPLATE(T, mat_t) X,
                const TEMPLATE(T, mat_t) A, const TEMPLATE(T, mat_t) B,
                                                 const TEMPLATE(T, ctx_t) ctx);

/* Transforms ****************************************************************/


void TEMPLATE(T, mat_similarity) (TEMPLATE(T, mat_t) A, slong r,
                               TEMPLATE(T, t) d, const TEMPLATE(T, ctx_t) ctx);

/* Characteristic polynomial *************************************************/

/* this prototype really lives in fq_poly_templates.h
 * FQ_MAT_TEMPLATES_INLINE
 * void TEMPLATE(T, mat_charpoly)(TEMPLATE(T, poly_t) p,
 *                          TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx)
 * {
 *   TEMPLATE(T, mat_charpoly_danilevsky) (p, A, ctx);
 * }
 */

/* Minimal polynomial ********************************************************/

/* this prototype really lives in fq_poly_templates.h
 *
 * void TEMPLATE(T, mat_minpoly) (TEMPLATE(T, poly_t) p,
 *                   const TEMPLATE(T, mat_t) X, const TEMPLATE(T, ctx_t) ctx);
 */

#ifdef __cplusplus
}
#endif

#endif

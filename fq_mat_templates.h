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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#ifdef T

#include "flint.h"
#include "templates.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    TEMPLATE(T, struct) * entries;
    slong r;
    slong c;
    TEMPLATE(T, struct) ** rows;
} TEMPLATE(T, mat_struct);

typedef TEMPLATE(T, mat_struct) TEMPLATE(T, mat_t)[1];

/* Memory management  ********************************************************/

FLINT_DLL void TEMPLATE(T, mat_init)(TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_init_set)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, mat_t) src,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_swap)(TEMPLATE(T, mat_t) mat1, TEMPLATE(T, mat_t) mat2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_set)(TEMPLATE(T, mat_t) mat1, const TEMPLATE(T, mat_t) mat2,
                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_clear)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, mat_equal)(const TEMPLATE(T, mat_t) mat1,
                       const TEMPLATE(T, mat_t) mat2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, mat_is_zero)(const TEMPLATE(T, mat_t) mat,
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

FQ_MAT_TEMPLATES_INLINE TEMPLATE(T, struct) *
TEMPLATE(T, mat_entry)(const TEMPLATE(T, mat_t) mat, slong i, slong j)
{
    return mat->rows[i] + j;
}

FQ_MAT_TEMPLATES_INLINE void
TEMPLATE(T, mat_entry_set)(TEMPLATE(T, mat_t) mat, slong i, slong j,
                           const TEMPLATE(T, t) x,
                           const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, set)(TEMPLATE(T, mat_entry)(mat, i, j), x, ctx);
}

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

/* Assignment  ***************************************************************/

FLINT_DLL void TEMPLATE(T, mat_zero)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

/* Windows and concatenation */

FLINT_DLL void TEMPLATE(T, mat_window_init)(TEMPLATE(T, mat_t) window,
                             const TEMPLATE(T, mat_t) mat,
                             slong r1, slong c1, slong r2, slong c2,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_window_clear)(TEMPLATE(T, mat_t) window,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_concat_horizontal)(TEMPLATE(T, mat_t) res,
                           const TEMPLATE(T, mat_t) mat1,  const TEMPLATE(T, mat_t) mat2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_concat_vertical)(TEMPLATE(T, mat_t) res,
                           const TEMPLATE(T, mat_t) mat1,  const TEMPLATE(T, mat_t) mat2,
                           const TEMPLATE(T, ctx_t) ctx);


/* Input and output  *********************************************************/

FLINT_DLL int TEMPLATE(T, mat_fprint)(FILE * file, const TEMPLATE(T, mat_t) mat,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, mat_fprint_pretty)(FILE * file, const TEMPLATE(T, mat_t) mat,
                                   const TEMPLATE(T, ctx_t) ctx);

FQ_MAT_TEMPLATES_INLINE
int TEMPLATE(T, mat_print)(const TEMPLATE(T, mat_t) mat,
                           const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, mat_fprint)(stdout, mat, ctx);
}

FQ_MAT_TEMPLATES_INLINE
int TEMPLATE(T, mat_print_pretty)(const TEMPLATE(T, mat_t) mat,
                                  const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, mat_fprint_pretty)(stdout, mat, ctx);
}

/* TODO: Read functions */

/* Random matrix generation  *************************************************/

FLINT_DLL void TEMPLATE(T, mat_randtest)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_randrank)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          slong rank, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, mat_randpermdiag)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                              TEMPLATE(T, struct) * diag, slong n,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_randops)(TEMPLATE(T, mat_t) mat, slong count,
                         flint_rand_t state, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_randtril)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          int unit, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_randtriu)(TEMPLATE(T, mat_t) mat, flint_rand_t state,
                          int unit, const TEMPLATE(T, ctx_t) ctx);

/* Norms */

/* Transpose */

/* Addition and subtraction */

FLINT_DLL void TEMPLATE(T, mat_add)(TEMPLATE(T, mat_t) C,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_sub)(TEMPLATE(T, mat_t) C,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_neg)(TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, mat_t) A,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_submul)(TEMPLATE(T, mat_t) D,
                        const TEMPLATE(T, mat_t) C,
                        const TEMPLATE(T, mat_t) A,
                        const TEMPLATE(T, mat_t) B,
                        const TEMPLATE(T, ctx_t) ctx);

/* Scalar operations */

/* Multiplication */

FLINT_DLL void TEMPLATE(T, mat_mul)(TEMPLATE(T, mat_t) C,
                     const TEMPLATE(T, mat_t) A,
                     const TEMPLATE(T, mat_t) B,
                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_mul_classical)(TEMPLATE(T, mat_t) C,
                               const TEMPLATE(T, mat_t) A,
                               const TEMPLATE(T, mat_t) B,
                               const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL void TEMPLATE(T, mat_mul_KS)(TEMPLATE(T, mat_t) C,
                        const TEMPLATE(T, mat_t) A,
                        const TEMPLATE(T, mat_t) B,
                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, mat_lu)(slong * P,
                    TEMPLATE(T, mat_t) A,
                    int rank_check,
                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, mat_lu_recursive)(slong * P,
                              TEMPLATE(T, mat_t) A,
                              int rank_check,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, mat_lu_classical)(slong * P, TEMPLATE(T, mat_t) A, int rank_check,
                              const TEMPLATE(T, ctx_t) ctx);


/* Solving *******************************************************************/

FLINT_DLL slong TEMPLATE(T, mat_rref)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, mat_nullspace)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) A,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, mat_rank)(const TEMPLATE(T, mat_t) A,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_solve_tril)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) L,
                            const TEMPLATE(T, mat_t) B, int unit,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_solve_tril_classical)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) L,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_solve_tril_recursive)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) L,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_solve_triu)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) U,
                            const TEMPLATE(T, mat_t) B, int unit,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_solve_triu_classical)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) U,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, mat_solve_triu_recursive)(TEMPLATE(T, mat_t) X,
                                      const TEMPLATE(T, mat_t) U,
                                      const TEMPLATE(T, mat_t) B,
                                      int unit,
                                      const TEMPLATE(T, ctx_t) ctx);



#ifdef __cplusplus
}
#endif

#endif

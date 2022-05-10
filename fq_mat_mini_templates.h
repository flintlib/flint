/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL void TEMPLATE(T, mat_init)(TEMPLATE(T, mat_t) mat, slong rows, slong cols,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_clear)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, ctx_t) ctx);

FQ_MAT_MINI_TEMPLATES_INLINE TEMPLATE(T, struct) *
TEMPLATE(T, mat_entry)(const TEMPLATE(T, mat_t) mat, slong i, slong j)
{
    return mat->rows[i] + j;
}

/* assignment  ****************************************************************/

FLINT_DLL void TEMPLATE(T, mat_set)(TEMPLATE(T, mat_t) mat1, const TEMPLATE(T, mat_t) mat2,
                     const TEMPLATE(T, ctx_t) ctx);

FQ_MAT_MINI_TEMPLATES_INLINE
void
TEMPLATE(T, mat_init_set)(TEMPLATE(T, mat_t) mat, const TEMPLATE(T, mat_t) src,
                          const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_init) (mat, src->r, src->c, ctx);
    TEMPLATE(T, mat_set) (mat, src, ctx);
}

FLINT_DLL void TEMPLATE(T, mat_zero)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mat_one)(TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

/* arithmetic operations ******************************************************/

FLINT_DLL void TEMPLATE(T, mat_mul)(TEMPLATE(T, mat_t) C,
                     const TEMPLATE(T, mat_t) A,
                     const TEMPLATE(T, mat_t) B,
                     const TEMPLATE(T, ctx_t) ctx);

/* solving ********************************************************************/

FLINT_DLL int TEMPLATE(T, mat_solve)(TEMPLATE(T, mat_t) X,
               const TEMPLATE(T, mat_t A), const TEMPLATE(T, mat_t) C,
                                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, mat_can_solve)(TEMPLATE(T, mat_t) X,
                const TEMPLATE(T, mat_t) A, const TEMPLATE(T, mat_t) B,
                                                 const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif

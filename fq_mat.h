/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_MAT_H
#define FQ_MAT_H

#ifdef FQ_MAT_INLINES_C
#define FQ_MAT_TEMPLATES_INLINE FLINT_DLL
#define FQ_MAT_INLINE FLINT_DLL
#else
#define FQ_MAT_TEMPLATES_INLINE static __inline__
#define FQ_MAT_INLINE static __inline__
#endif


#include "fq.h"
#include "fq_vec.h"

/* Cutoff between classical and recursive triangular solving */
#define FQ_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define FQ_MAT_SOLVE_TRI_COLS_CUTOFF 64

/* Cutoff between classical and recursive LU decomposition */
#define FQ_MAT_LU_RECURSIVE_CUTOFF 4

FQ_MAT_INLINE int FQ_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_ctx_t ctx)
{
    if (5 * FLINT_MIN(r, c) > 8 * fq_ctx_degree(ctx) + 29)
        return 1;
    else
        return 0;
}

#define T fq
#define CAP_T FQ
#include "fq_mat_templates.h"
#undef CAP_T
#undef T

#endif

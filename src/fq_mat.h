/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_MAT_H
#define FQ_MAT_H

#ifdef FQ_MAT_INLINES_C
#define FQ_MAT_TEMPLATES_INLINE
#define FQ_MAT_INLINE
#else
#define FQ_MAT_TEMPLATES_INLINE static inline
#define FQ_MAT_INLINE static inline
#endif

#include "fq_types.h"

/* Cutoff between classical and recursive triangular solving */
#define FQ_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define FQ_MAT_SOLVE_TRI_COLS_CUTOFF 64

/* Cutoff between classical and recursive LU decomposition */
#define FQ_MAT_LU_RECURSIVE_CUTOFF 4

int FQ_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_ctx_t ctx);

#define T fq
#define CAP_T FQ
#include "fq_mat_templates.h"
#undef CAP_T
#undef T

#endif

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

#ifndef FQ_NMOD_MAT_H
#define FQ_NMOD_MAT_H

#ifdef FQ_NMOD_MAT_INLINES_C
#define FQ_MAT_TEMPLATES_INLINE FLINT_DLL
#define FQ_NMOD_MAT_INLINE FLINT_DLL
#else
#define FQ_MAT_TEMPLATES_INLINE static __inline__
#define FQ_NMOD_MAT_INLINE static __inline__
#endif

#include "fq_nmod.h"
#include "fq_nmod_vec.h"

/* Cutoff between classical and recursive triangular solving */
#define FQ_NMOD_MAT_SOLVE_TRI_ROWS_CUTOFF 64
#define FQ_NMOD_MAT_SOLVE_TRI_COLS_CUTOFF 64

/* Cutoff between classical and recursive LU decomposition */
#define FQ_NMOD_MAT_LU_RECURSIVE_CUTOFF 4

FQ_NMOD_MAT_INLINE
int FQ_NMOD_MAT_MUL_KS_CUTOFF(slong r, slong c, const fq_nmod_ctx_t ctx)
{
    if (FLINT_MIN(r, c) > fq_nmod_ctx_degree(ctx) / 20 + 6)
        return 1;
    else
        return 0;
}

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_mat_templates.h"
#undef CAP_T
#undef T

#endif

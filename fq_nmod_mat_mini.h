/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_MAT_MINI_H
#define FQ_NMOD_MAT_MINI_H

#ifdef FQ_NMOD_MAT_INLINES_C
#define FQ_MAT_TEMPLATES_INLINE FLINT_DLL
#define FQ_NMOD_MAT_INLINE FLINT_DLL
#else
#define FQ_MAT_TEMPLATES_INLINE static __inline__
#define FQ_NMOD_MAT_INLINE static __inline__
#endif

#include "flint.h"

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_mat_mini_templates.h"
#undef CAP_T
#undef T

#endif

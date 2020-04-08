/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FQ_SPARSE_VEC_H
#define FQ_SPARSE_VEC_H

#ifdef FQ_SPARSE_VEC_INLINES_C
#define FQ_SPARSE_VEC_TEMPLATES_INLINE FLINT_DLL
#define FQ_SPARSE_VEC_INLINE FLINT_DLL
#else
#define FQ_SPARSE_VEC_TEMPLATES_INLINE static __inline__
#define FQ_SPARSE_VEC_INLINE static __inline__
#endif


#include "fq.h"
#include "fq_vec.h"

#define T fq
#define CAP_T FQ
#include "fq_sparse_vec_templates.h"
#undef CAP_T
#undef T

#endif

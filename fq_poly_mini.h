/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_POLY_MINI_H
#define FQ_POLY_MINI_H

#ifdef FQ_POLY_INLINES_C
#define FQ_POLY_INLINE FLINT_DLL
#define FQ_POLY_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_POLY_INLINE static __inline__
#define FQ_POLY_TEMPLATES_INLINE static __inline__
#endif

#include "fq_mini.h"

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_poly_mini_templates.h"
#undef CAP_T
#undef T

#endif

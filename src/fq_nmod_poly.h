/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_POLY_H
#define FQ_NMOD_POLY_H

#ifdef FQ_NMOD_POLY_INLINES_C
#define FQ_POLY_TEMPLATES_INLINE
#define FQ_NMOD_POLY_INLINE
#else
#define FQ_POLY_TEMPLATES_INLINE static inline
#define FQ_NMOD_POLY_INLINE static inline
#endif

#include "fq_nmod_types.h"

#define FQ_NMOD_POLY_DIVREM_DIVCONQUER_CUTOFF  16
#define FQ_NMOD_COMPOSE_MOD_LENH_CUTOFF 6
#define FQ_NMOD_COMPOSE_MOD_PREINV_LENH_CUTOFF 6
#define FQ_NMOD_MUL_CLASSICAL_CUTOFF 6
#define FQ_NMOD_SQR_CLASSICAL_CUTOFF 6
#define FQ_NMOD_MULLOW_CLASSICAL_CUTOFF 6

#define FQ_NMOD_POLY_HGCD_CUTOFF 25
#define FQ_NMOD_POLY_SMALL_GCD_CUTOFF 110
#define FQ_NMOD_POLY_GCD_CUTOFF 120

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_poly_templates.h"
#undef CAP_T
#undef T

#endif

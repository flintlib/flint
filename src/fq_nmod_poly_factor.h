/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_POLY_FACTOR_H
#define FQ_NMOD_POLY_FACTOR_H

#ifdef FQ_NMOD_POLY_FACTOR_INLINES_C
#define FQ_POLY_FACTOR_TEMPLATES_INLINE
#define FQ_NMOD_POLY_FACTOR_INLINE
#else
#define FQ_POLY_FACTOR_TEMPLATES_INLINE static inline
#define FQ_NMOD_POLY_FACTOR_INLINE static inline
#endif

#include "fq_nmod_types.h"

int FQ_NMOD_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_nmod_ctx_t ctx, slong length);

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#include "fq_poly_factor_templates.h"
#undef CAP_T
#undef T

/* Inlines *******************************************************************/

void fq_nmod_poly_factor_get_poly(fq_nmod_poly_t z,
            const fq_nmod_poly_factor_t fac, slong i, const fq_nmod_ctx_t ctx);

#endif

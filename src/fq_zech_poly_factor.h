/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_POLY_FACTOR_H
#define FQ_ZECH_POLY_FACTOR_H

#ifdef FQ_ZECH_POLY_FACTOR_INLINES_C
#define FQ_POLY_FACTOR_TEMPLATES_INLINE
#define FQ_ZECH_POLY_FACTOR_INLINE
#else
#define FQ_POLY_FACTOR_TEMPLATES_INLINE static inline
#define FQ_ZECH_POLY_FACTOR_INLINE static inline
#endif

#include "fq_zech_types.h"

int FQ_ZECH_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_zech_ctx_t ctx, slong length);

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#include "fq_poly_factor_templates.h"
#undef CAP_T
#undef T

void fq_zech_poly_factor_get_poly(fq_zech_poly_t z,
            const fq_zech_poly_factor_t fac, slong i, const fq_zech_ctx_t ctx);

#endif

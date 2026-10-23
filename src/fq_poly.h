/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_POLY_H
#define FQ_POLY_H

#ifdef FQ_POLY_INLINES_C
#define FQ_POLY_INLINE
#define FQ_POLY_TEMPLATES_INLINE
#else
#define FQ_POLY_INLINE static inline
#define FQ_POLY_TEMPLATES_INLINE static inline
#endif

#include "fq_types.h"

#define FQ_POLY_DIVREM_DIVCONQUER_CUTOFF  16
#define FQ_COMPOSE_MOD_LENH_CUTOFF 6
#define FQ_COMPOSE_MOD_PREINV_LENH_CUTOFF 6
#define FQ_MUL_CLASSICAL_CUTOFF 6
#define FQ_MULLOW_CLASSICAL_CUTOFF 6
#define FQ_SQR_CLASSICAL_CUTOFF 6

#define FQ_POLY_HGCD_CUTOFF 30
#define FQ_POLY_SMALL_GCD_CUTOFF 80
#define FQ_POLY_GCD_CUTOFF 90

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#include "fq_poly_templates.h"
#undef CAP_T
#undef T

/* not available for fq_nmod and fq_zech */
#ifdef __cplusplus
extern "C" {
#endif

void _fq_poly_mul_reorder(fq_struct * rop, const fq_struct * op1, slong len1,
                          const fq_struct * op2, slong len2, const fq_ctx_t ctx);
void fq_poly_mul_reorder(fq_poly_t rop, const fq_poly_t op1,
                         const fq_poly_t op2, const fq_ctx_t ctx);

void _fq_poly_sqr_reorder(fq_struct * rop, const fq_struct * op, slong len,
                          const fq_ctx_t ctx);
void fq_poly_sqr_reorder(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

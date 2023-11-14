/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_EMBED_H
#define FQ_EMBED_H

#ifdef FQ_EMBED_INLINES_C
#define FQ_EMBED_INLINE
#define FQ_EMBED_TEMPLATES_INLINE
#else
#define FQ_EMBED_INLINE static inline
#define FQ_EMBED_TEMPLATES_INLINE static inline
#endif

#include "fq_types.h"

#define T fq
#define B fmpz_mod
#include "fq_embed_templates.h"
#undef B
#undef T

void fq_modulus_pow_series_inv(fmpz_mod_poly_t res, const fq_ctx_t ctx, slong trunc);

void fq_modulus_derivative_inv(fq_t m_prime, fq_t m_prime_inv, const fq_ctx_t ctx);

#endif

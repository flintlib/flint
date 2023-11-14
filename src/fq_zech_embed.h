/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_EMBED_H
#define FQ_ZECH_EMBED_H

#ifdef FQ_ZECH_EMBED_INLINES_C
#define FQ_ZECH_EMBED_INLINE
#define FQ_EMBED_TEMPLATES_INLINE
#else
#define FQ_ZECH_EMBED_INLINE static inline
#define FQ_EMBED_TEMPLATES_INLINE static inline
#endif

#include "fq_zech_types.h"

#define T fq_zech
#define B nmod
#include "fq_embed_templates.h"

void TEMPLATE(T, modulus_pow_series_inv)(TEMPLATE(B, poly_t) res,
                                         const TEMPLATE(T, ctx_t) ctx,
                                         slong trunc);

#undef B
#undef T

void fq_zech_modulus_derivative_inv(fq_zech_t m_prime, fq_zech_t m_prime_inv,
                                                     const fq_zech_ctx_t ctx);

#endif

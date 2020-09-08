/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_NMOD_EMBED_H
#define FQ_NMOD_EMBED_H

#ifdef FQ_NMOD_EMBED_INLINES_C
#define FQ_NMOD_EMBED_INLINE FLINT_DLL
#define FQ_EMBED_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_NMOD_EMBED_INLINE static __inline__
#define FQ_EMBED_TEMPLATES_INLINE static __inline__
#endif

#include "fq_nmod.h"

#define T fq_nmod
#define B nmod
#include "fq_embed_templates.h"

FQ_EMBED_TEMPLATES_INLINE
void TEMPLATE(T, modulus_pow_series_inv)(TEMPLATE(B, poly_t) res,
                                         const TEMPLATE(T, ctx_t) ctx,
                                         slong trunc)
{
    TEMPLATE(B, poly_reverse)(res, 
                              TEMPLATE(T, ctx_modulus)(ctx), 
                              TEMPLATE(T, ctx_degree)(ctx) + 1);
    TEMPLATE(B, poly_inv_series)(res, res, trunc);
}

#undef B
#undef T

FQ_NMOD_EMBED_INLINE void fq_nmod_modulus_derivative_inv(fq_nmod_t m_prime,
                                                         fq_nmod_t m_prime_inv,
                                                         const fq_nmod_ctx_t ctx)
{
    nmod_poly_derivative(m_prime, fq_nmod_ctx_modulus(ctx));
    fq_nmod_inv(m_prime_inv, m_prime, ctx);
}

#endif

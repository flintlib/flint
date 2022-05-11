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
#define FQ_EMBED_INLINE FLINT_DLL
#define FQ_EMBED_TEMPLATES_INLINE FLINT_DLL
#else
#define FQ_EMBED_INLINE static __inline__
#define FQ_EMBED_TEMPLATES_INLINE static __inline__
#endif

#include "fq_mini.h"
#include "fmpz_mod_poly_mini.h"

/* Define some macros here in order to avoid inclusions of other headers. */
#define FQ_CTX_MODULUS(ctx)     \
    ((ctx)->modulus)

#define FQ_CTX_DEGREE(ctx)      \
    ((ctx)->modulus->length - 1)

#define T fq
#define CAP_T FQ
#define B fmpz_mod
#include "fq_embed_templates.h"
#undef B
#undef CAP_T
#undef T

FQ_EMBED_INLINE
void fq_modulus_pow_series_inv(fmpz_mod_poly_t res, const fq_ctx_t ctx, slong trunc)
{
    fmpz_mod_poly_reverse(res, FQ_CTX_MODULUS(ctx), FQ_CTX_DEGREE(ctx) + 1, ctx->ctxp);
    fmpz_mod_poly_inv_series(res, res, trunc, ctx->ctxp);
}

FQ_EMBED_INLINE
void fq_modulus_derivative_inv(fq_t m_prime, fq_t m_prime_inv, const fq_ctx_t ctx)
{
    fmpz_mod_poly_t tmp;
    fmpz_mod_poly_init(tmp, ctx->ctxp);
    
    fmpz_mod_poly_derivative(tmp, FQ_CTX_MODULUS(ctx), ctx->ctxp);
    fmpz_poly_fit_length(m_prime, tmp->length);
    _fmpz_vec_set(m_prime->coeffs, tmp->coeffs, tmp->length);
    _fmpz_poly_set_length(m_prime, tmp->length);
    
    fq_inv(m_prime_inv, m_prime, ctx);
    fmpz_mod_poly_clear(tmp, ctx->ctxp);
}

#endif

/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fq.h"
#include "fq_embed.h"

void fq_modulus_derivative_inv(fq_t m_prime, fq_t m_prime_inv, const fq_ctx_t ctx)
{
    fmpz_mod_poly_t tmp;
    fmpz_mod_poly_init(tmp, ctx->ctxp);

    fmpz_mod_poly_derivative(tmp, fq_ctx_modulus(ctx), ctx->ctxp);
    fmpz_poly_fit_length(m_prime, tmp->length);
    _fmpz_vec_set(m_prime->coeffs, tmp->coeffs, tmp->length);
    _fmpz_poly_set_length(m_prime, tmp->length);

    fq_inv(m_prime_inv, m_prime, ctx);
    fmpz_mod_poly_clear(tmp, ctx->ctxp);
}

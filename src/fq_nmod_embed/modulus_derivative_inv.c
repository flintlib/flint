/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "fq_nmod.h"
#include "fq_nmod_embed.h"

void fq_nmod_modulus_derivative_inv(fq_nmod_t m_prime, fq_nmod_t m_prime_inv,
                                                     const fq_nmod_ctx_t ctx)
{
    nmod_poly_derivative(m_prime, fq_nmod_ctx_modulus(ctx));
    fq_nmod_inv(m_prime_inv, m_prime, ctx);
}

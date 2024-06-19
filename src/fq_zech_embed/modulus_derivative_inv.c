/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "fq_nmod_embed.h"
#include "fq_zech.h"
#include "fq_zech_embed.h"

void fq_zech_modulus_derivative_inv(fq_zech_t m_prime, fq_zech_t m_prime_inv,
                                                     const fq_zech_ctx_t ctx)
{
    fq_nmod_t m_nmod, m_inv_nmod;
    fq_nmod_init(m_nmod, ctx->fq_nmod_ctx);
    fq_nmod_init(m_inv_nmod, ctx->fq_nmod_ctx);

    fq_nmod_modulus_derivative_inv(m_nmod, m_inv_nmod, ctx->fq_nmod_ctx);

    fq_zech_set_fq_nmod(m_prime, m_nmod, ctx);
    fq_zech_set_fq_nmod(m_prime_inv, m_inv_nmod, ctx);

    fq_nmod_clear(m_nmod, ctx->fq_nmod_ctx);
    fq_nmod_clear(m_inv_nmod, ctx->fq_nmod_ctx);
}

/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly_factor.h"


void _fq_nmod_mpoly_get_lead0(
    fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpolyl_lead_coeff(c, A, 1, ctx);
}

void _fq_nmod_mpoly_set_lead0(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong deg;
    fq_nmod_mpoly_t t, g;

    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_init(g, ctx);

    deg = fq_nmod_mpoly_degree_si(B, 0, ctx);
    FLINT_ASSERT(deg >= 0);
    fq_nmod_mpoly_gen(g, 0, ctx);
    fq_nmod_mpoly_pow_ui(g, g, deg, ctx);
    _fq_nmod_mpoly_get_lead0(t, B, ctx);
    fq_nmod_mpoly_sub(t, c, t, ctx);
    fq_nmod_mpoly_mul(t, t, g, ctx);
    fq_nmod_mpoly_add(A, B, t, ctx);

    fq_nmod_mpoly_clear(t, ctx);
    fq_nmod_mpoly_clear(g, ctx);
}

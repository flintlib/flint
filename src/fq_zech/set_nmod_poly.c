/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void fq_zech_set_nmod_poly(fq_zech_t a, const nmod_poly_t b, const fq_zech_ctx_t ctx)
{
    ulong blen = b->length;
    const mp_limb_t * bcoeffs = b->coeffs;
    mp_limb_t qm1 = ctx->qm1;
    ulong i;
    fq_zech_t t;
    fq_zech_zero(a, ctx);

    for (i = 0; i < blen; i++)
    {
        if (bcoeffs[i] == 0)
            continue;
        t->value = (blen <= qm1) ? i : (i % qm1);
        fq_zech_mul_ui(t, t, bcoeffs[i], ctx);
        fq_zech_add(a, a, t, ctx);
    }
}

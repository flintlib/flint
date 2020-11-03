/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_fq_poly_neg(
    n_poly_t A,
    const n_poly_t B,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Blen = B->length;

    n_poly_fit_length(A, d*Blen);
    _nmod_vec_neg(A->coeffs, B->coeffs, d*Blen, ctx->mod);
    A->length = Blen;
    _n_fq_poly_normalise(A, d);
}

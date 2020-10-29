/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_fq_poly_add(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Blen = B->length;
    slong Clen = C->length;

    if (Blen > Clen)
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        if (A != B)
            _nmod_vec_set(A->coeffs + d*Clen, B->coeffs + d*Clen, d*(Blen - Clen));
        A->length = Blen;
    }
    else if (Blen < Clen)
    {
        n_poly_fit_length(A, d*Clen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Blen, ctx->mod);
        if (A != C)
            _nmod_vec_set(A->coeffs + d*Blen, C->coeffs + d*Blen, d*(Clen - Blen));
        A->length = Clen;
    }
    else
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_add(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        A->length = Clen;
        _n_fq_poly_normalise(A, d);
    }

    FLINT_ASSERT(n_fq_poly_is_canonical(A, ctx));
}

/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_fq_poly_sub(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    const n_fq_poly_t C,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    slong Blen = B->length;
    slong Clen = C->length;

    if (Blen > Clen)
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_sub(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        if (A != B)
            for (i = d*Clen; i < d*Blen; i++)
                A->coeffs[i] = B->coeffs[i];
        A->length = Blen;
    }
    else if (Blen < Clen)
    {
        n_poly_fit_length(A, d*Clen);
        _nmod_vec_sub(A->coeffs, B->coeffs, C->coeffs, d*Blen, ctx->mod);
        for (i = d*Blen; i < d*Clen; i++)
            A->coeffs[i] = nmod_neg(C->coeffs[i], ctx->mod);
        A->length = Clen;
    }
    else
    {
        n_poly_fit_length(A, d*Blen);
        _nmod_vec_sub(A->coeffs, B->coeffs, C->coeffs, d*Clen, ctx->mod);
        A->length = Clen;
        _n_fq_poly_normalise(A, d);
    }
}

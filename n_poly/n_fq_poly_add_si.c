/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"


void n_fq_poly_add_si(
    n_fq_poly_t A,
    const n_fq_poly_t B,
    slong c,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);

    if (A != B)
        n_fq_poly_set(A, B, ctx);

    if (A->length < 1)
    {
        n_poly_fit_length(A, d);
        A->length = 1;
    }

    n_fq_add_si(A->coeffs + d*0, A->coeffs + d*0, c, ctx);

    _n_fq_poly_normalise(A, d);
}

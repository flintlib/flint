/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"

void fq_nmod_get_nmod_mat(nmod_mat_t col,
                          const fq_nmod_t a,
                          const fq_nmod_ctx_t ctx)
{
    slong i, n = fq_nmod_ctx_degree(ctx);
    for (i = 0; i < a->length; i++)
        nmod_mat_set_entry(col, i, 0, a->coeffs[i]);
    for ( ; i < n; i++)
        nmod_mat_entry(col, i, 0) = 0;
}

void fq_nmod_set_nmod_mat(fq_nmod_t a,
                          const nmod_mat_t col,
                          const fq_nmod_ctx_t ctx)
{
    slong i, n = fq_nmod_ctx_degree(ctx);
    nmod_poly_fit_length(a, n);
    a->length = n;
    for (i = 0; i < n; i++)
        a->coeffs[i] = nmod_mat_entry(col, i, 0);
    _nmod_poly_normalise(a);
}

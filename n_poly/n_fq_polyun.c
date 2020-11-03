/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"

void n_fq_polyun_set(n_fq_polyun_t A, const n_fq_polyun_t B, const fq_nmod_ctx_t ctx)
{
    slong i;
    n_polyun_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        A->terms[i].exp = B->terms[i].exp;
        n_fq_poly_set(A->terms[i].coeff, B->terms[i].coeff, ctx);
    }
    A->length = B->length;
}


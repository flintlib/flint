/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void fq_get_fmpz_mod_mat(fmpz_mod_mat_t col,
                         const fq_t a,
                         const fq_ctx_t ctx)
{
    slong i, n = fq_ctx_degree(ctx);
    for (i = 0; i < a->length; i++)
        fmpz_mod_mat_set_entry(col, i, 0, a->coeffs + i);
    for ( ; i < n; i++)
        fmpz_zero(fmpz_mod_mat_entry(col, i, 0));
}

void fq_set_fmpz_mod_mat(fq_t a,
                         const fmpz_mod_mat_t col,
                         const fq_ctx_t ctx)
{
    slong i, n = fq_ctx_degree(ctx);
    fmpz_poly_fit_length(a, n);
    a->length = n;
    for (i = 0; i < n; i++)
        fmpz_set(a->coeffs + i, fmpz_mod_mat_entry(col, i, 0));
    _fmpz_poly_normalise(a);
}

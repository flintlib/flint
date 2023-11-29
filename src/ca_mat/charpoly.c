/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
_ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    if (ca_mat_nrows(mat) <= 2)
    {
        _ca_mat_charpoly_berkowitz(cp, mat, ctx);
    }
    else
    {
/*
        ca_field_ptr K;

        K = _ca_mat_same_field(mat, ctx);

        if (0 && K != NULL && CA_FIELD_IS_NF(K))
        {
            if (_ca_mat_charpoly_danilevsky(cp, mat, ctx))
                return;
        }
*/

        _ca_mat_charpoly_berkowitz(cp, mat, ctx);
    }
}

void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)
{
    if (mat->r != mat->c)
    {
        flint_throw(FLINT_ERROR, "Exception (ca_mat_charpoly).  Non-square matrix.\n");
    }

    ca_poly_fit_length(cp, mat->r + 1, ctx);
    _ca_poly_set_length(cp, mat->r + 1, ctx);
    _ca_mat_charpoly(cp->coeffs, mat, ctx);
}

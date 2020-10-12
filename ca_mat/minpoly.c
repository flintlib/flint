/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int ca_mat_minpoly(ca_poly_t res, const ca_mat_t mat, ca_ctx_t ctx)
{
    if (mat->r != mat->c)
    {
        flint_printf("Exception (ca_mat_minpoly).  Non-square matrix.\n");
        flint_abort();
    }

    ca_mat_charpoly(res, mat, ctx);
    return ca_poly_squarefree_part(res, res, ctx);
}

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
ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx)
{
    if (mat->entries != NULL)
    {
        _ca_vec_clear(mat->entries, mat->r * mat->c, ctx);
        flint_free(mat->rows);
    }
}

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
ca_mat_window_init(ca_mat_t window, const ca_mat_t mat,
    slong r1, slong c1, slong r2, slong c2, ca_ctx_t ctx)
{
    slong i;
    window->entries = NULL;

    window->rows = flint_malloc((r2 - r1) * sizeof(ca_ptr));

    for (i = 0; i < r2 - r1; i++)
        window->rows[i] = mat->rows[r1 + i] + c1;

    window->r = r2 - r1;
    window->c = c2 - c1;
}

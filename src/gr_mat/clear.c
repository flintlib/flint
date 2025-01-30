/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

void
gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx)
{
    if (mat->entries != NULL)
    {
        /* Should always be the case for non-window matrices. */
        FLINT_ASSERT(mat->stride == mat->c);

        _gr_vec_clear(mat->entries, mat->r * mat->c, ctx);
        flint_free(mat->entries);
    }
}

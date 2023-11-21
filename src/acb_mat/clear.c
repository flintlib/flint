/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_clear(acb_mat_t mat)
{
    if (mat->entries != NULL)
    {
        _acb_vec_clear(mat->entries, mat->r * mat->c);
        flint_free(mat->rows);
    }
}

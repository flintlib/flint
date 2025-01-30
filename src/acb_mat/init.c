/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"

void
acb_mat_init(acb_mat_t mat, slong r, slong c)
{
    mat->entries = NULL;
    mat->r = r;
    mat->c = c;
    mat->stride = c;

    if (r != 0 && c != 0)
        mat->entries = _acb_vec_init(r * c);
}

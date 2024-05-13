/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_mat.h"

ulong
nmod_mat_trace(const nmod_mat_t mat)
{
    ulong t;
    slong i, n = nmod_mat_nrows(mat);

    if (n == 0)
        return 0;

    t = nmod_mat_entry(mat, 0, 0);

    for (i = 1; i < n; i++)
        t = nmod_add(t, nmod_mat_entry(mat, i, i), mat->mod);

    return t;
}

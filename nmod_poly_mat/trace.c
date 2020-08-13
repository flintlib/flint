/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_trace(nmod_poly_t trace, const nmod_poly_mat_t mat)
{
    slong i, n = nmod_poly_mat_nrows(mat);

    if (n == 0)
        nmod_poly_zero(trace);
    else
    {
        nmod_poly_set(trace, nmod_poly_mat_entry(mat, 0, 0));
        for (i = 1; i < n; i++)
            nmod_poly_add(trace, trace, nmod_poly_mat_entry(mat, i, i));
    }
}

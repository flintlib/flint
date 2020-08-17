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
nmod_poly_mat_one(nmod_poly_mat_t A)
{
    slong i, n;

    nmod_poly_mat_zero(A);
    n = FLINT_MIN(A->r, A->c);

    for (i = 0; i < n; i++)
        nmod_poly_one(nmod_poly_mat_entry(A, i, i));
}

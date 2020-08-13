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
nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A)
{
    slong n = A->r;

    if (n == 0)
    {
        nmod_poly_one(det);
    }
    else if (n == 1)
    {
        nmod_poly_set(det, nmod_poly_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        nmod_poly_t tmp;
        nmod_poly_init(tmp, nmod_poly_mat_modulus(A));
        nmod_poly_mul(det, nmod_poly_mat_entry(A, 0, 0),
                           nmod_poly_mat_entry(A, 1, 1));
        nmod_poly_mul(tmp, nmod_poly_mat_entry(A, 0, 1),
                           nmod_poly_mat_entry(A, 1, 0));
        nmod_poly_sub(det, det, tmp);
        nmod_poly_clear(tmp);
    }
    else if (n < 15)  /* should be entry sensitive too */
    {
        nmod_poly_mat_det_fflu(det, A);
    }
    else
    {
        nmod_poly_mat_det_interpolate(det, A);
    }
}

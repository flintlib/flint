/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_det(nmod_poly_t det, const nmod_poly_mat_t A)
{
    len_t n = A->r;

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

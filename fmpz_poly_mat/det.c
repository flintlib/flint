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

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_det(fmpz_poly_t det, const fmpz_poly_mat_t A)
{
    len_t n = A->r;

    if (n == 0)
    {
        fmpz_poly_set_ui(det, 1UL);
    }
    else if (n == 1)
    {
        fmpz_poly_set(det, fmpz_poly_mat_entry(A, 0, 0));
    }
    else if (n == 2)
    {
        fmpz_poly_t tmp;
        fmpz_poly_init(tmp);
        fmpz_poly_mul(det, fmpz_poly_mat_entry(A, 0, 0),
                           fmpz_poly_mat_entry(A, 1, 1));
        fmpz_poly_mul(tmp, fmpz_poly_mat_entry(A, 0, 1),
                           fmpz_poly_mat_entry(A, 1, 0));
        fmpz_poly_sub(det, det, tmp);
        fmpz_poly_clear(tmp);
    }
    else if (n < 15)  /* should be entry sensitive too */
    {
        fmpz_poly_mat_det_fflu(det, A);
    }
    else
    {
        fmpz_poly_mat_det_interpolate(det, A);
    }
}

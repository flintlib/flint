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
fmpz_poly_mat_det_interpolate(fmpz_poly_t det, const fmpz_poly_mat_t A)
{
    len_t i, l, n, len;

    fmpz_mat_t X;
    fmpz * x;
    fmpz * d;

    n = A->r;

    if (n == 0)
    {
        fmpz_poly_one(det);
        return;
    }

    l = fmpz_poly_mat_max_length(A);

    if (l == 0)
    {
        fmpz_poly_zero(det);
        return;
    }

    /* Bound degree based on Laplace expansion */
    len = n*(l - 1) + 1;

    x = _fmpz_vec_init(len);
    d = _fmpz_vec_init(len);
    fmpz_mat_init(X, n, n);

    for (i = 0; i < len; i++)
    {
        fmpz_set_si(x + i, -len/2 + i);
        fmpz_poly_mat_evaluate_fmpz(X, A, x + i);
        fmpz_mat_det(d + i, X);
    }

    fmpz_poly_interpolate_fmpz_vec(det, x, d, len);

    _fmpz_vec_clear(x, len);
    _fmpz_vec_clear(d, len);
    fmpz_mat_clear(X);
}

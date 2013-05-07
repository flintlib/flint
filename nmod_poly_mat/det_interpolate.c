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

#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_det_interpolate(nmod_poly_t det, const nmod_poly_mat_t A)
{
    long i, l, n, len;

    nmod_mat_t X;
    mp_ptr x, d;

    n = A->r;

    if (n == 0)
    {
        nmod_poly_one(det);
        return;
    }

    l = nmod_poly_mat_max_length(A);

    if (l == 0)
    {
        nmod_poly_zero(det);
        return;
    }

    /* Bound degree based on Laplace expansion */
    len = n*(l - 1) + 1;

    /* Not enough points to interpolate */
    if (len > nmod_poly_mat_modulus(A))
    {
        nmod_poly_mat_det_fflu(det, A);
        return;
    }

    x = _nmod_vec_init(len);
    d = _nmod_vec_init(len);
    nmod_mat_init(X, n, n, nmod_poly_mat_modulus(A));

    for (i = 0; i < len; i++)
    {
        x[i] = i;
        nmod_poly_mat_evaluate_nmod(X, A, x[i]);
        d[i] = nmod_mat_det(X);
    }

    nmod_poly_interpolate_nmod_vec(det, x, d, len);

    _nmod_vec_clear(x);
    _nmod_vec_clear(d);
    nmod_mat_clear(X);
}

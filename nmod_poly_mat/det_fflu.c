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
#include "perm.h"

void
nmod_poly_mat_det_fflu(nmod_poly_t det, const nmod_poly_mat_t A)
{
    len_t n = nmod_poly_mat_nrows(A);

    if (n == 0)
        nmod_poly_one(det);
    else
    {
        nmod_poly_mat_t tmp;
        len_t * perm;
        nmod_poly_mat_init_set(tmp, A);
        perm = _perm_init(n);

        nmod_poly_mat_fflu(tmp, det, perm, tmp, 1);
        if (_perm_parity(perm, n))
            nmod_poly_neg(det, det);

        _perm_clear(perm);
        nmod_poly_mat_clear(tmp);
    }
}

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
#include "perm.h"

void
fmpz_poly_mat_det_fflu(fmpz_poly_t det, const fmpz_poly_mat_t A)
{
    len_t n = fmpz_poly_mat_nrows(A);

    if (n == 0)
        fmpz_poly_one(det);
    else
    {
        fmpz_poly_mat_t tmp;
        len_t * perm;
        fmpz_poly_mat_init_set(tmp, A);
        perm = _perm_init(n);

        fmpz_poly_mat_fflu(tmp, det, perm, tmp, 1);
        if (_perm_parity(perm, n))
            fmpz_poly_neg(det, det);

        _perm_clear(perm);
        fmpz_poly_mat_clear(tmp);
    }
}

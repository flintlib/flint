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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"

void fmpz_sparse_set_coeff(fmpz_sparse_t poly, const fmpz_t c, const fmpz_t e)
{
    slong ind;
    ind = _fmpz_sparse_index(poly, e);
    if (ind < 0) {
        if (! fmpz_is_zero(c)) {
            /* nonzero coeff, not in the polynomial; we have to insert it. */
            _fmpz_sparse_reserve(poly, poly->length+1);
            ind = -1 - ind;
            _fmpz_sparse_vec_shift(poly, ind, poly->length, 1);
            fmpz_set(poly->coeffs+ind, c);
            fmpz_set(poly->expons+ind, e);
            ++ poly->length;
        }
    }
    else if (fmpz_is_zero(c)) {
        /* zero coeff, in the polynomial; we have to remove it. */
        fmpz_clear(poly->coeffs+ind);
        fmpz_clear(poly->expons+ind);
        _fmpz_sparse_vec_shift(poly, ind+1, poly->length, -1);
        poly->length -= 1;
    }
    else {
        /* nonzero coeff, in the polynomial; just change it. */
        fmpz_set(poly->coeffs + ind, c);
    }
}

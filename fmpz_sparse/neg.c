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

void fmpz_sparse_neg(fmpz_sparse_t res, const fmpz_sparse_t poly)
{
    if (res == poly) {
        int i;
        for (i=0; i<poly->length; ++i) fmpz_neg(poly->coeffs + i, poly->coeffs + i);
    }
    else {
        /* similar to set(), but you negate. */
        int cur=0;
        for (; cur < res->length && cur < poly->length; ++cur) {
            fmpz_neg(res->coeffs+cur, poly->coeffs+cur);
            fmpz_neg(res->expons+cur, poly->expons+cur);
        }
        for (; cur < res->length; ++cur) {
            fmpz_clear(res->coeffs+cur);
            fmpz_clear(res->expons+cur);
        }
        _fmpz_sparse_reserve(res, poly->length);
        res->length = poly->length;
        for (; cur < poly->length; ++cur) {
            fmpz_init(res->coeffs+cur);
            fmpz_init(res->expons+cur);
            fmpz_neg(res->coeffs+cur, poly->coeffs+cur);
            fmpz_neg(res->expons+cur, poly->expons+cur);
        }
    }
}

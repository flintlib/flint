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

void fmpz_sparse_set(fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    int cur=0;
    for (; cur < poly1->length && cur < poly2->length; ++cur) {
        fmpz_set(poly1->coeffs+cur, poly2->coeffs+cur);
        fmpz_set(poly1->expons+cur, poly2->expons+cur);
    }
    for (; cur < poly1->length; ++cur) {
        fmpz_clear(poly1->coeffs+cur);
        fmpz_clear(poly1->expons+cur);
    }
    _fmpz_sparse_reserve(poly1, poly2->length);
    poly1->length = poly2->length;
    for (; cur < poly2->length; ++cur) {
        fmpz_init(poly1->coeffs+cur);
        fmpz_init(poly1->expons+cur);
        fmpz_set(poly1->coeffs+cur, poly2->coeffs+cur);
        fmpz_set(poly1->expons+cur, poly2->expons+cur);
    }
}

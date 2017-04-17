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

void _fmpz_sparse_append_si(fmpz_sparse_t poly1, const fmpz_sparse_t poly2, 
    slong c, ulong e)
{
    slong curf, curg;
    _fmpz_sparse_reserve(poly1, poly1->length + poly2->length);
    for (curf = poly1->length, curg=0; curg < poly2->length; ++curf, ++curg) {
        fmpz_init(poly1->coeffs + curf);
        fmpz_init(poly1->expons + curf);
        fmpz_mul_si(poly1->coeffs + curf, poly2->coeffs + curg, c);
        fmpz_add_ui(poly1->expons + curf, poly2->expons + curg, e);
    }
    poly1->length += poly2->length;
}

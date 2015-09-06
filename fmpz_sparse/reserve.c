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

void _fmpz_sparse_reserve(fmpz_sparse_t poly, slong terms)
{
    if (terms > poly->alloc) {
        if (terms < 2*poly->alloc) terms = 2*poly->alloc;
        poly->coeffs = (fmpz*) flint_realloc(poly->coeffs, terms*sizeof(fmpz));
        poly->expons = (fmpz*) flint_realloc(poly->expons, terms*sizeof(fmpz));
        poly->alloc = terms;
    }
}

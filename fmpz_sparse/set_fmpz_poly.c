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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_poly.h"

void fmpz_sparse_set_fmpz_poly(fmpz_sparse_t poly1, const fmpz_poly_t poly2)
{
    slong i, j = 0;
    fmpz_sparse_zero(poly1);
    _fmpz_sparse_reserve(poly1, poly2->length);
    for (i = fmpz_poly_degree(poly2); i>=0; --i) 
    {
      if (!fmpz_is_zero(poly2->coeffs + i)) 
      {
        fmpz_init(poly1->expons + j);
        fmpz_init(poly1->coeffs + j);
        fmpz_set_si(poly1->expons + j, i + 1);
        fmpz_set(poly1->coeffs + j, poly2->coeffs + i);
        j += 1;
      }
    }
    poly1->length = j;
    _fmpz_sparse_normalise(poly1);
}

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

void fmpz_sparse_set_fmpz_poly(fmpz_sparse_t poly1, const fmpz_poly_t poly2)
{
    slong i;
    fmpz_t c, e;
    fmpz_init(c);
    fmpz_init(e);
    fmpz_sparse_zero(poly1);
    for (i = fmpz_poly_degree(poly2); i>=0; --i) {
        fmpz_poly_get_coeff_fmpz(c, poly2, i);
        if (!fmpz_is_zero(c)) {
            fmpz_set_si(e, i);
            fmpz_sparse_set_coeff(poly1, c, e);
        }
    }
}

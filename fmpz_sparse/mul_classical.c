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

void fmpz_sparse_mul_classical(fmpz_sparse_t res, 
    const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
    slong curg;
    if (res == poly1 || res == poly2) {
        fmpz_sparse_t temp;
        fmpz_sparse_init(temp);
        fmpz_sparse_mul_classical(temp, poly1, poly2);
        fmpz_sparse_set(res, temp);
        fmpz_sparse_clear(temp);
        return;
    }
    /* TODO rewrite this if addition is rewritten. */
    fmpz_sparse_zero(res);
    for (curg=0; curg < poly2->length; ++curg) {
        _fmpz_sparse_append(res, poly1, poly2->coeffs+curg, poly2->expons+curg);
    }
    _fmpz_sparse_normalise(res);
}

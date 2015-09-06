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

void fmpz_sparse_rem_cyc_dense(fmpz_poly_t res, const fmpz_sparse_t poly, ulong e)
{
    slong i;
    fmpz_poly_zero(res);
    fmpz_poly_set_coeff_si(res, e, 1); /* set to x^e, to reserve zeros */
    for (i=0; i < poly->length; ++i) {
        slong rese = fmpz_fdiv_ui(poly->expons + i, e);
        fmpz* resc = fmpz_poly_get_coeff_ptr(res, rese);
        fmpz_add(resc, resc, poly->coeffs + i);
    }
    fmpz_poly_set_coeff_si(res, e, 0); /* set x^e coeff back to 0 */
}

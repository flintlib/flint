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

void fmpz_sparse_get_fmpz_poly(fmpz_poly_t out, const fmpz_sparse_t in)
{
    slong i;
    fmpz_t c, e;
    fmpz_init(c);
    fmpz_init(e);
    fmpz_poly_zero(out);
    for (i=0; i<fmpz_sparse_terms(in); ++i) {
        fmpz_sparse_get_term(c, e, in, i);
        FLINT_ASSERT (fmpz_fits_si(e));
        fmpz_poly_set_coeff_fmpz(out, fmpz_get_si(e), c);
    }
}

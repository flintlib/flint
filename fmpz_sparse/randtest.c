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

void fmpz_sparse_randtest(fmpz_sparse_t res, flint_rand_t state, 
    slong terms, const fmpz_t degree, mp_bitcnt_t bits)
{
    slong i;
    fmpz_sparse_zero(res);
    
    if(fmpz_is_zero(degree))
      return;

    if(fmpz_equal_si(degree, terms))
    {
      terms = fmpz_get_si(degree);
    }

    _fmpz_sparse_reserve(res, terms);
    for (i=0; i<terms; ++i) {
        fmpz_init(res->coeffs+i);
        fmpz_init(res->expons+i);
        fmpz_randbits(res->coeffs+i, state, bits);
        fmpz_randm(res->expons+i, state, degree);
    }
    res->length = terms;
    _fmpz_sparse_normalise(res);
}
